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

function slicetimings(img, lags, mon_cyc, nslices::Int)
    @assert ndims(img) == 3 # should be stored by Imagine as a 2D timeseries
    #first set of images is always the target / template image
    target = img[:,:,1:nslices]
    testimgs = view(img, :, :, (nslices+1):size(img,3))
    fwdimgs = view(testimgs, :, :, 1:2:size(testimgs,3))
    backimgs = view(testimgs, :, :, 2:2:size(testimgs,3))
    nlags =length(lags)
    ntrials = convert(Int, size(fwdimgs,3) / (nlags*nslices))
    fwd_timings = eltype(lags)[]
    back_timings = eltype(lags)[]
    @showprogress for i = 1:nslices
        fixed = view(target, :, :, i)
        fwd_lag, fwd_mm = best_lag(fixed, view(fwdimgs, :, :, (i-1)*(ntrials*nlags)+1:i*(ntrials*nlags)), lags)
        push!(fwd_timings, fwd_lag)
        back_lag, back_mm = best_lag(fixed, view(backimgs, :, :, (i-1)*(ntrials*nlags)+1:i*(ntrials*nlags)), lags)
        push!(back_timings, back_lag)
        print("Slice $i complete\n")
        @show fwd_mm
        @show back_mm
    end
    return fwd_timings, back_timings
end

#best_lag may be a misnomer here because the lag we choose doesn't correspond with the sensor lag and it's not observed to be
#consistent per-slice.  It depends on so many factors that we just do it emprically.
function best_lag(target, testimgs, lags)
    nlags = length(lags)
    ntrials = convert(Int, size(testimgs,3) / nlags)
    bestlag = lags[1]
    best_mm = Inf
	maxradians = pi/10
	rgridsz = 7
	mxshift = (20,20)
    alg = RigidGridStart(copy(target), maxradians, rgridsz, mxshift; print_level=5, max_iter=200)
    for i = 1:length(lags)
        moving = squeeze(mean(view(testimgs, :, :, (i-1)*ntrials+1:i*ntrials), 3),3)
		mon = monitor(alg, ())
		mon[:tform] = nothing
		mon[:mismatch] = 0.0
		mon = driver(alg, moving, mon)
		mm = mon[:mismatch]
		if mm < best_mm
			best_mm = mm
			bestlag = lags[i]
		end
	end
	@show best_mm
	@show bestlag
	return bestlag, best_mm
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
