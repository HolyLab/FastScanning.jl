#Convenience for extracting positioner mod and mon with error checks
function pos_mod_mon(coms, recs)
    pos = getpositioners(coms)
    if length(pos) > 1 error("Provide only one positioner command.") end
    pos_mon = getpositionermonitors(recs)
    if length(pos_mon) > 1 error("Provide only one positioner monitor.") end
    return pos[1], pos_mon[1]
end

#take the mean MON cycle and calculate the number of sampls that MON is shifted from MOD
function mon_lag_nsamps(pos, pos_mon, ncycs_ignore)
    pos_cyc_name = first(sequence_names(pos))
    pos_cycle = get_samples(pos, pos_cyc_name)
    @assert all(sequence_names(pos).==pos_cyc_name)
    nsamps_cycle = Base.length(pos_cycle) #should always be even
    mon_samps = get_samples(pos_mon)
    mon_samps = view(mon_samps, (ncycs_ignore*nsamps_cycle+1):Base.length(mon_samps))
    cycs_mon = ImagineAnalyses.get_cycles(mon_samps, nsamps_cycle)
    mean_cyc = mean(ustrip.(cycs_mon), 2)[:] * unit(mon_samps[1])
    return ImagineAnalyses.mon_delay(pos_cycle, mean_cyc)
end

function local_minima(mms)
    min_is = Int[]
    if mms[2] > mms[1]
        push!(mms, 1)
    end
    for i = 2:(length(mms)-1)
        if mms[i-1] > mms[i] < mms[i+1]
            push!(min_is, i)
        end
    end
    if mms[end] < mms[end-1]
        push!(min_is,length(mms))
        warn("Found local minimum at the last temporal offset.  You should use a wider range of offsets")
    end
    return min_is
end

#Instead of just taking the lowest mismatch, find all local minima and take the highest-indexed minimum that satisfies thresh_fac
function ind_best_offset(mms; thresh_fac = 0.1)
    min_is = local_minima(mms)
    min_mm = minimum(mms)
    #thresh = min_mm + (maximum(mms) - min_mm) * thresh_fac
    thresh = min_mm + min_mm * thresh_fac
    for i = length(min_is):-1:1
        if mms[min_is[i]] <= thresh
            return min_is[i]
        end
    end
end

#imageseq is a 3D image array where each 2D slice is an image acquired for calibration
#imageseq should be either forward or reverse slice trials, but not both (use exclude_fwd to separate)
#The first (slow) stack should not be included in imgseq (use exclude_template_stack to get that one)
get_slice_trials(imgseq, sliceidx, ntrials, ntoffsets) = view(imgseq, :, :, (sliceidx-1)*(ntrials*ntoffsets)+1:sliceidx*(ntrials*ntoffsets))
exclude_template_stack(imgseq, nslices) = view(imgseq, :, :, (nslices+1):size(imgseq,3))
get_template_stack(imgseq, nslices) = imgseq[:,:,1:nslices]
#assumes template stack is excluded
get_fwd_imgs(imgseq) = view(imgseq, :, :, 1:2:size(imgseq,3))
get_bck_imgs(imgseq) = view(imgseq, :, :, 2:2:size(imgseq,3))

#NOTE: the timings returned here are additional offsets relative to any lag already applied with a Calibration when 
#generating the command signals.

function slicetimings(img, toffsets, nslices::Int; allow_shifts=true, allow_rotations=false)
    timings, mms, tfms = timings_mms_tfms(img, toffsets, nslices; allow_shifts=allow_shifts, allow_rotations=allow_rotations)
    return timings
end

function timings_mms_tfms(img, toffsets, nslices; allow_shifts=true, allow_rotations=false, record_all=false)
    @assert ndims(img) == 3 # should be stored by Imagine as a 2D timeseries
    #first set of images is always the target / template image
    target = get_template_stack(img, nslices)
    testimgs = exclude_template_stack(img, nslices)
    fwdimgs = get_fwd_imgs(testimgs)
    backimgs = get_bck_imgs(testimgs)
    timings_mms_tfms(target, fwdimgs, backimgs, toffsets, nslices; allow_shifts=allow_shifts, allow_rotations=allow_rotations, record_all=record_all)
end

#Setting record_all will return arrays-of-arrays of mismatches and transforms (one for each condition) rather than just the best
function timings_mms_tfms(target, fwdimgs, backimgs, toffsets, nslices; allow_shifts=true, allow_rotations=false, record_all=false)
    ntoffsets =length(toffsets)
    ntrials = convert(Int, size(fwdimgs,3) / (ntoffsets*nslices))
    fwd_timings = eltype(toffsets)[]
    back_timings = eltype(toffsets)[]
    fwd_tfms = []; back_tfms = []; fwd_mms = []; back_mms = []
    @showprogress for i = 1:nslices
        @show i
        fixed = view(target, :, :, i)
        movfwd = get_slice_trials(fwdimgs, i, ntrials, ntoffsets)
        ftfms, fmms = toffset_fits(fixed, movfwd, toffsets; allow_shifts=allow_shifts, allow_rotations=allow_rotations)
        #fwdi = ind_best_offset(fmms)
        fwdi = indmin(fmms)
        @show fwd_toffset = toffsets[fwdi]
        @show fwd_mm = fmms[fwdi]
        push!(fwd_timings, fwd_toffset)
        movback = get_slice_trials(backimgs, i, ntrials, ntoffsets)
        btfms, bmms = toffset_fits(fixed, movback, toffsets; allow_shifts=allow_shifts, allow_rotations = allow_rotations)
        #backi = ind_best_offset(bmms)
        backi = indmin(bmms)
        @show back_toffset = toffsets[backi]
        @show back_mm = bmms[backi]
        push!(back_timings, back_toffset)
        if record_all
            push!(fwd_tfms, ftfms)
            push!(back_tfms, btfms)
            push!(fwd_mms, fmms)
            push!(back_mms, bmms)
        else
            push!(fwd_tfms, ftfms[fwdi])
            push!(back_tfms, btfms[backi])
            push!(fwd_mms, fwd_mm)
            push!(back_mms, back_mm)
        end
        print("Slice $i complete\n")
    end
    return (fwd_timings, back_timings), (fwd_mms, back_mms), (fwd_tfms, back_tfms)
end

function preprocess(img::AbstractArray{Float64,2}; sigmas=(1.0,1.0), sqrt_tfm=true, bias=Float64(reinterpret(N0f16, UInt16(100))))
    img = imfilter(img, KernelFactors.IIRGaussian(Float64, sigmas))
    output = similar(img)
    for i in eachindex(img)
        temp = max(0.0, img[i] - bias)
        output[i] = ifelse(sqrt_tfm, sqrt(temp), temp)
    end
    return output
end

function align2d(fixed, moving; thresh_fac=0.9, sigmas=(1.0,1.0), allow_shifts = true, allow_rotations =false)
    mxshift = (8,8)
    moving = preprocess(moving; sigmas=sigmas, sqrt_tfm=false)
    if allow_rotations
        maxradians = pi/180
        thresh = (1-thresh_fac) * sum(abs2.(fixed[.!(isnan.(fixed))]))
        tfm, mm = qd_rigid(fixed, moving, mxshift, [maxradians], [pi/10000]; thresh=thresh, tfm0=IdentityTransformation())
        return tfm, mm
#        rgridsz = 7
#        alg = RigidGridStart(fixed, maxradians, rgridsz, mxshift; thresh_fac=thresh_fac, print_level=0, max_iter=100)
#        mon = monitor(alg, ())
#        mon[:tform] = nothing
#        mon[:mismatch] = 0.0
#        mon = driver(alg, moving, mon)
#        return mon[:tform], mon[:mismatch]
    elseif allow_shifts
        thresh = (1-thresh_fac) * sum(abs2.(fixed[.!(isnan.(fixed))]))
        shft, mm = RegisterOptimize.best_shift(fixed, moving, mxshift, thresh; normalization=:intensity, initial_tfm=IdentityTransformation())
        return tformtranslate([shft...]), mm
    else
        mm = mismatch0(fixed, moving; normalization=:intensity)
        return tformtranslate([0;0]), ratio(mm, 0.0)
    end
end

#get_offset_trials returns the set of trials at a particular temporal offset
#Assumes that img_seq:
#   -does not include the slow "target" stack
#   -includes only forward or reverse images, not both
#   -includes data for only one slice
#(This means img_seq has size ntrials * ntoffsets in dimension 3)
get_toffset_trials(img_seq, itoffset, ntrials) = view(img_seq, :, :, (itoffset-1)*ntrials+1:itoffset*ntrials)

#consistent per-slice.  It depends on so many factors that we just do it emprically.
function toffset_fits(target, testimgs, toffsets::Vector; sigmas=(1.0,1.0), allow_shifts=true, allow_rotations = false)
    ntoffsets = length(toffsets)::Int
    ntrials = convert(Int, size(testimgs,3) / ntoffsets)
    tfms = []
    mms = Float64[]
    target = preprocess(Float64.(target); sigmas=sigmas, sqrt_tfm=false)
    #Threads.@threads for i = 1:length(toffsets)  #Currently segfaults on julia 0.6.2
    for i = 1:length(toffsets)
        moving = squeeze(mean(Float64.(get_toffset_trials(testimgs, i, ntrials)), 3),3)
        tfm, mm = align2d(target, moving; sigmas=sigmas, allow_shifts=allow_shifts, allow_rotations=allow_rotations)
        #@show mm
        push!(tfms, tfm)
        push!(mms, mm)
    end
    return tfms, mms
end

function find_bidi_centers(cyc::AbstractVector, sampr::HasInverseTimeUnits)
    @show pmin = minimum(cyc)
    @show pmax = maximum(cyc)
    pctr = pmin + (pmax-pmin)/2 #we will place the flash in the center
    find_bidi_locations(cyc, sampr, pctr)
end

function find_bidi_locations(cyc::AbstractVector, sampr::HasInverseTimeUnits, ploc::HasLengthUnits)
    stack_rate = 2 * (sampr / length(cyc))
    pad_secs = (1/stack_rate)/200 #wait 0.5% of a cycle before accepting threshold crossings (helps with noise)
    pad_nsamps = round(Int, pad_secs * sampr)
    mon_is = ImagineAnalyses.find_circular(cyc, [ploc;], pad_nsamps)[1]
    @assert Base.length(mon_is) == 2 #should have just two crossings, one fwd one back
    fwd_i = mon_is[1]
    back_i = mon_is[2]
    return fwd_i, back_i
end
