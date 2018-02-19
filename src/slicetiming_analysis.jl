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

#NOTE: the timings returned here are additional offsets relative to any lag already applied with a Calibration when 
#generating the command signals.

function slicetimings(img, toffsets, mon_cyc, nslices::Int; allow_rotations=false)
    @assert ndims(img) == 3 # should be stored by Imagine as a 2D timeseries
    #first set of images is always the target / template image
    target = img[:,:,1:nslices]
    testimgs = view(img, :, :, (nslices+1):size(img,3))
    fwdimgs = view(testimgs, :, :, 1:2:size(testimgs,3))
    backimgs = view(testimgs, :, :, 2:2:size(testimgs,3))
    ntoffsets =length(toffsets)
    ntrials = convert(Int, size(fwdimgs,3) / (ntoffsets*nslices))
    fwd_timings = eltype(toffsets)[]
    back_timings = eltype(toffsets)[]
    @showprogress for i = 1:nslices
        @show i
        fixed = view(target, :, :, i)
        movfwd = view(fwdimgs, :, :, (i-1)*(ntrials*ntoffsets)+1:i*(ntrials*ntoffsets))
        ftfms, fmms = toffset_fits(fixed, movfwd, toffsets; allow_rotations = allow_rotations)
        @show fwd_toffset = toffsets[indmin(fmms)]
        @show fwd_mm = minimum(fmms)
        push!(fwd_timings, fwd_toffset)
        movback = view(backimgs, :, :, (i-1)*(ntrials*ntoffsets)+1:i*(ntrials*ntoffsets))
        btfms, bmms = toffset_fits(fixed, movback, toffsets; allow_rotations = allow_rotations)
        @show back_toffset = toffsets[indmin(bmms)]
        @show back_mm = minimum(bmms)
        push!(back_timings, back_toffset)
        print("Slice $i complete\n")
    end
    return fwd_timings, back_timings
end

function align2d(fixed, moving; thresh_fac=0.9, sigmas=(1.0,1.0), allow_rotations =false)
    mxshift = (16,16)
    moving = imfilter(Float64.(moving), KernelFactors.IIRGaussian(Float64, sigmas))
    fixed = imfilter(Float64.(fixed), KernelFactors.IIRGaussian(Float64, sigmas))
    if allow_rotations
        maxradians = pi/180
        rgridsz = 7
        alg = RigidGridStart(fixed, maxradians, rgridsz, mxshift; thresh_fac=thresh_fac, print_level=0, max_iter=100)
        mon = monitor(alg, ())
        mon[:tform] = nothing
        mon[:mismatch] = 0.0
        mon = driver(alg, moving, mon)
        return mon[:tform], mon[:mismatch]
    else
        thresh = (1-thresh_fac) * sum(abs2.(fixed[.!(isnan.(fixed))]))
        shft, mm = RegisterOptimize.best_shift(fixed, moving, mxshift, thresh; normalization=:intensity, initial_tfm=IdentityTransformation())
        return tformtranslate([shft...]), mm
	end
end

#consistent per-slice.  It depends on so many factors that we just do it emprically.
function toffset_fits(target, testimgs, toffsets; sigmas=(1.0,1.0), allow_rotations = false)
    ntoffsets = length(toffsets)
    ntrials = convert(Int, size(testimgs,3) / ntoffsets)
    tfms = AffineTransform[]
    mms = Float64[]
    for i = 1:length(toffsets)
        moving = squeeze(mean(Float64.(view(testimgs, :, :, (i-1)*ntrials+1:i*ntrials)), 3),3)
        tfm, mm = align2d(target, moving; sigmas=sigmas, allow_rotations=allow_rotations)
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
