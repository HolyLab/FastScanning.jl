#using Unitful, ImagineInterface, ImagineAnalyses, Images, ProgressMeter
#import Unitful:s, μm, V
#
##Generate an ouput ImagineSignal for commanding the piezo to repeat the sample sequence `one_cyc` `ncycles` times.
##Appends those samples to `pos` and also returns dummy camera and laser signals so that Imagine doesn't complain
##`ncycles` default to 20 about seconds-worth of cycles based on the sample rate of `pos`
#function pos_experiment(pos::ImagineSignal, one_cyc::Vector, ncycles=ceil(Int, 20.0 / ustrip(inv(samprate(pos))*length(one_cyc))))
#    @assert isempty(pos)
#    rig_sigs = rigtemplate(rig_name(pos); sample_rate = samprate(pos))
#    add_sequence!(pos, "bidi_cycle", pos)
#    cycle_nsamps = length(one_cyc)
#    cam = first(getcameras(rig_sigs)) #not used
#    las = first(getlasers(rig_sigs)) #not used
#    add_sequence!(cam, "quiet", falses(cycle_nsamps)) #just to pass validation
#    for i = 1:ncycles
#        append!(pos, "bidi_cycle")
#        append!(cam, "quiet")
#        append!(las, "quiet")
#    end
#    return [pos;cam;las]
#end
#
#function pos_experiment(rig::AbstractString, pos_name::AbstractString, pstart, pstop, stack_rate, ncycles=ceil(Int, 20.0/ustrip(inv(stack_rate))); sample_rate = 100000Hz)
#    rig_sigs = rigtemplate(rig; sample_rate = sample_rate)
#    posfwd, posback = gen_bidi_pos(pstart, pstop, 1/stack_rate, sample_rate)
#    pos = getname(rig_sigs, pos_name)
#    pos_experiment(pos, vcat(posfwd, posback), sample_rate, ncycles)
#end
#
##take the mean MON cycle and calculate the number of sampls that MON is shifted from MOD
#function mon_lag_nsamps(pos, pos_mon, ncycs_ignore)
#    pos_cyc_name = first(sequence_names(pos))
#    pos_cycle = get_samples(pos, pos_cyc_name)
#    @assert all(sequence_names(pos).==pos_cyc_name)
#    nsamps_cycle = Base.length(pos_cycle) #should always be even
#    mon_samps = get_samples(pos_mon)
#    mon_samps = view(mon_samps, (ncycs_ignore*nsamps_cycle+1):Base.length(mon_samps))
#    cycs_mon = ImagineAnalyses.get_cycles(mon_samps, nsamps_cycle)
#    mean_cyc = mean(ustrip.(cycs_mon), 2)[:] * unit(mon_samps[1])
#    return ImagineAnalyses.mon_delay(pos_cycle, mean_cyc)
#end
#
##Convenience for extracting positioner mod and mon with error checks
#function pos_mod_mon(coms, recs)
#    pos = getpositioners(coms)
#    if length(pos) > 1 error("Provide only one positioner command.") end
#    pos_mon = getpositionermonitors(recs)
#    if length(pos_mon) > 1 error("Provide only one positioner monitor.") end
#    return pos[1], pos_mon[1]
#end
#
##reads piezo monitor samples from recs and returns an imagine procedure with the following attributes:
##The SignalGenerator generates command signals with two phases:
##   Phase 1:  A slow (0.2Hz) stack is acquired with camera and laser pulses timed to match target z locations determined by z_spacing.
##   Phase 2:  Many piezo cycles with two camera and laser pulses per cycle.  The piezo cycle is taken from the recs input.
##           Each pulse is timed at a different lag relative to the target z location determined by z_spacing (and used to place slices in Phase 1)
##The process function expects the recorded signals from the above commands as inputs, and it outputs time indices indicating when to time slices so that z location
##when imaging fast matches the location of images from the template stack of phase 1 above.
##Optionally pass a Calibration to be applied to the piezo MON signal before doing all this (if, for example, you know the constant lag in the MON output)
#function slicetiming_experiment(coms, recs, las_name::AbstractString, cam_name::AbstractString, z_spacing::HasLengthUnits, z_pad::HasLengthUnits; cal = Calibration())
#    pos, pos_mon = pos_mod_mon(coms, recs)
#    pos_mon_cal = apply(cal, pos_mon)
#    mon_samps = get_samples(pos_mon_cal)
#    nsamps_cycle = full_length(first(sequences(pos)))
#    ncycs_ignore = ceil(Int, (ustrip(samprate(pos))*5) / nsamps_cycle) #5 seconds of warmup
#    mean_cyc = mean_cycle(mon_samps, nsamps_cycle; start_idx = ncycs_ignore*nsamps_cycle+1)
#    #set pmin and pmax empirically
#    pmin = minimum(mean_cyc)
#    pmax = maximum(mean_cyc)
#    slice_zs = ImagineInterface.slice_positions(pmin, pmax, z_spacing, z_pad)
#    slicetiming_experiment(pos, mean_cyc, las_name, cam_name, slice_zs)
#end
#
#function slicetiming_experiment(pos::ImagineSignal, mon_cyc::ImagineSignal, las_name::AbstractString, cam_name::AbstractString, slice_zs)
#    #@show lags = [0.0002s:0.00005s:0.0007s...]  #These worked to find the lag on the analog piezo of OCPI2
#    @show lags = [0.0005s:0.00005s:0.002s...] #These worked to find the lag on the digital piezo of OCPI2
#    gen_desc = "Generate a set camera and laser commands to find slice timings empirically given a completed dynamic positioner recording"
#    sig_gen = SignalGenerator(gen_desc, slicetiming_experiment, (pos, mon_cyc, las_name, cam_name, lags, slice_zs))
#    desc = "Returns slice timings to align forward and reverse stacks at the desired spacing for high-speed volume acquisitions"
#    analyze_f = (sigs, img) -> slicetimings(sigs, img, lags_pre, mon_cyc, length(slice_zs))
#    #analyze_f = (coms, img) -> calc_lag(coms, img, lags_pre, ncycs_ignore, mon_cyc, nsamps_offset)
#    return ImagineProcedure(desc, sig_gen, analyze_f)
#end
#
##Similar to slicetiming_experiment but only checks one slice per stack, and uses the result to estimate the lag between the piezo's true position and the MON signal
##function pos_mon_lag_experiment(pos::ImagineSignal, mon_cyc::ImagineSignal, las_name::AbstractString, cam_name::AbstractString, slice_zs)
#
#function pos_mon_lag_experiment(pos::ImagineSignal, mon_cyc, las_name::AbstractString, cam_name::AbstractString)
#    #set pmin and pmax empirically
#    pmin = minimum(mon_cyc)
#    pmax = maximum(mon_cyc)
#    pctr = [(pmax-pmin)/2 * pmin]
#    #@show lags = [0.0002s:0.00005s:0.0007s...]  #These worked to find the lag on the analog piezo of OCPI2
#    @show lags = [0.0005s:0.00005s:0.002s...] #These worked to find the lag on the digital piezo of OCPI2
#    gen_desc = "Generate a set camera and laser commands to find piezo monitor lag empirically given a completed dynamic positioner recording"
#    sig_gen = SignalGenerator(gen_desc, slicetiming_experiment, (pos, mon_cyc, las_name, cam_name, lags, [pctr]))
#    desc = "Returns mean of two measurements: the monitor lag when taking an image while sweeping forward and also while sweeping back"
#    analyze_f = (sigs, img) -> mean(slicetimings(sigs, img, lags_pre, mon_cyc, 1)) #take the mean of the forward and reverse lag (hopefully they are about equal)
#    return ImagineProcedure(desc, sig_gen, analyze_f)
#end
#
#pos_mon_lag_experiment(coms, recs, las_name::AbstractString, cam_name::AbstractString) = pos_mon_lag_experiment(pos_mod_mon(coms,recs)..., las_name, cam_name)
#
#function pos_mon_lag_experiment(pos::ImagineSignal, pos_mon::ImagineSignal, las_name::AbstractString, cam_name::AbstractString)
#    mon_samps = get_samples(pos_mon)
#    nsamps_cycle = full_length(first(sequences(pos)))
#    ncycs_ignore = ceil(Int, (ustrip(samprate(pos))*5) / nsamps_cycle) #5 seconds of warmup
#    mon_cyc = mean_cycle(mon_samps, nsamps_cycle; start_idx = ncycs_ignore*nsamps_cycle+1)
#    pos_mon_lag_experiment(pos, mon_cyc, las_name, cam_name)
#end
#
##NOTE: the timings returned here are additional offsets relative to any lag already applied with a Calibration when 
##generating the command signals.
#function slicetimings(sigs, img, lags, mon_cyc, nslices::Int)
#    @assert ndims(img) == 3 # should be stored by Imagine as a 2D timeseries
#    #first set of images is always the target / template image
#    target = img[:,:,1:nslices]
#    testimgs = view(img, :, :, (nslices+1):size(img,3))
#    fwdimgs = view(testimgs, :, :, 1:2:size(testimgs,3))
#    backimgs = view(testimgs, :, :, 2:2:size(testimgs,3))
#    nlags =length(lags)
#    ntrials = convert(Int, size(fwdimgs,3) / nlags)
#    fwd_timings = Int[]
#    back_timings = Int[]
#    @showprogress for i = 1:nslices
#        push!(fwd_timings, best_lag(view(target, :, :, i), view(fwdimgs, :, :, (i-1)*(ntrials*nlags)+1:i*(ntrials*nlags)), lags))
#        push!(back_timings, best_lag(view(target, :, :, i), view(backimgs, :, :, (i-1)*(ntrials*nlags)+1:i*(ntrials*nlags)), lags))
#    end
#    return vcat(fwd_timings, back_timings)
#end
#
##best_lag may be a misnomer here because the lag we choose doesn't correspond with the sensor lag and it's not observed to be
##consistent per-slice.  It depends on so many factors that we just do it emprically.
#function best_lag(target, testimgs, lags)
#    nlags = length(lags)
#    ntrials = convert(Int, size(testimgs,3) / nlags)
#    bestlag = lags[1]
#    best_mm = Inf
#	maxradians = pi/10
#	rgridsz = 7
#	mxshift = (20,20)
#	alg = RigidGridStart(target, maxradians, rgridsz, mxshift; print_level=5)
#    for i = 1:length(lags)
#        moving = squeeze(mean(view(testimgs, :, :, (i-1)*ntrials+1:i*ntrials), 3),3)
#		mon = monitor(alg, ())
#		mon[:tform] = nothing
#		mon[:mismatch] = 0.0
#		mon = driver(alg, moving, mon)
#		mm = mon[:mismatch]
#		if mm < best_mm
#			best_mm = mm
#			bestlag = lags[i]
#		end
#	end
#	@show best_mm
#	@show bestlag
#	return bestlag
#end
#
#function find_bidi_centers(cyc::AbstractVector, sampr::HasInverseTimeUnits)
#    @show pmin = minimum(cyc)
#    @show pmax = maximum(cyc)
#    pctr = pmin + (pmax-pmin)/2 #we will place the flash in the center
#	find_bidi_locations(cyc, sampr, pctr)
#end
#
#function find_bidi_locations(cyc::AbstractVector, sampr::HasInverseTimeUnits, ploc::HasLengthUnits)
#    stack_rate = 2 * (sampr / length(cyc))
#    pos_avg_speed = stack_rate * prange #in microns per second
#    pad_secs = (1/stack_rate)/2 #wait half a cycle before accepting threshold crossings (helps with noise)
#    pad_nsamps = round(Int, pad_secs * sampr)
#    mon_is = ImagineAnalyses.find_circular(cyc, [ploc;], pad_nsamps)[1]
#    @assert Base.length(mon_is) == 2 #should have just two crossings, one fwd one back
#    fwd_i = mon_is[1]
#    back_i = mon_is[2]
#    return fwd_i, back_i
#end
#
#function warmup!(pos, cam, las, high_las, pos_cycle, ncycs_warmup)
#    nsamps_cycle = Base.length(pos_cycle) #should always be even
#    add_sequence!(cam, "cam_warmup", falses(nsamps_cycle))
#    add_sequence!(las, "las_warmup", falses(nsamps_cycle))
#    add_sequence!(pos, "pos_warmup", pos_cycle)
#    if rig_name(pos) == "ocpi-2"
#        add_sequence!(high_las, "las_high_warmup", trues(nsamps_cycle))
#    end
#    for i = 1:ncycs_warmup
#        append!(pos, "pos_warmup")
#        append!(cam, "cam_warmup")
#        append!(las, "las_warmup")
#        append!(high_las, "las_high_warmup")
#    end
#	return pos, cam, las, high_las
#end
#
##delete below?
#    #After the ignored cycles each trial has this structure:
#    #pos: pos_cyc_name
#    #cam: "fwd_exp" -> "back_exp_$i" where i is an index into the lags vector
#    #las: "fwd_flash" -> "back_flash_$i"
#    #in order to know when/where to place the camera and laser pulses we must look at the pos_mon signal
#
#function get_trial(ctr_i, lag, sample_rate, nsamps_exp, nsamps_flash, nsamps_sweep, is_fwd::Bool)
#	nsamps_from_ctr = div(nsamps_flash,2)
#	nsamps_lag = ImagineInterface.calc_num_samps(lags, sample_rate)
#	if is_fwd
#		flash_itv = ClosedInterval(ctr_i-nsamps_lag-nsamps_from_ctr, ctr_i-nsamps_lag+nsamps_from_ctr)
#		exp_itv = ClosedInterval(ctr_i-nsamps_lag+nsamps_from_ctr-nsamps_exp-1, ctr_i-nsamps_lag+nsamps_from_ctr)
#	else
#		flash_itv = ClosedInterval(ctr_i-nsamps_lag-nsamps_from_ctr-nsamps_sweep, ctr_i-nsamps_lag+nsamps_from_ctr-nsamps_sweep)
#		exp_itv = ClosedInterval(ctr_i-nsamps_lag+nsamps_from_ctr-nsamps_exp-1-nsamps_sweep, ctr_i-nsamps_lag+nsamps_from_ctr-nsamps_sweep)
#	end
#	las_samps = ImagineInterface.gen_pulses(nsamps_sweep, [flash_itv;])
#	cam_samps = ImagineInterface.gen_pulses(nsamps_sweep, [exp_itv;])
#	return cam_samps, las_samps
#end
#
##returns the SignalGenerator for step2
#function slicetimings_step2_gen!(pos::T1, mod_cyc, mon_cyc, las, las_high, cam, lags, slice_zs) where {T1<:ImagineSignal}
#    nsamps_cycle = Base.length(mod_cyc) #should always be even
#    ncycs_ignore = ceil(Int, (ustrip(samprate(pos))*5) / nsamps_cycle) #5 seconds of warmup
#    @show nsamps_offset = mon_lag_nsamps(pos, mon_cyc, ncycs_ignore)
#	#temporal shift to make calculations easier (will shift things back later)
#    mon_cyc = circshift(mon_cyc, -nsamps_offset)
#    ntrials=10 #number of trials per lag value per z slice
#    #We want to deliver the laser pulse during the global shutter time (tglobal in PCO camera manual)
#    #If the user images with full chip then the global exposure segment doesn't begin until _10ms_ after the exposure starts.
#    #Therefore with a 1ms laser pulse we need an 11ms exposure time, and the pulse is placed in the last 1ms 
#    exp_time = 0.011s
#    flash_time = 0.001s
#    empty!(pos)
#	if nsamps_offset < 0
#		error("Cannot handle negative offset (but could with a bit more coding work)")
#	elseif nsamps_offset > 0
#		append!(las, "delay", fill(false, nsamps_offset))
#		append!(cam, "delay")
#		append!(high_las, "delayhigh", fill(true, nsamps_offset))
#	end
#    nsamps_flash = ImagineInterface.calc_num_samps(flash_time, samprate(pos))
#    if iseven(nsamps_flash)
#	    nsamps_flash -= 1 #force odd so that fwd and reverse flashes are easier to align
#    end
#	nsamps_exp = ImagineInterface.calc_num_samps(exp_time, samprate(pos))
#
#	#Sequence from here:
#	# one slow image stack
#	# several fast "warmup" piezo cycles with no camera / laser activity (ncycs_ignore)
#	# ntrials * nlags piezo cycles in which we take a pair of slices for each cycle
#	# repeat above step for each slice, so a total of ntrials * nlags * nslices cycles (and twice that many images)
#    ncycs_ignore = ceil(Int, (ustrip(samprate(pos))*5) / nsamps_cycle) #5 seconds of warmup
#	acquire_template!(pos, cam, las, high_las, slice_zs) #appends slow stack sequence to each signal input
#	warmup!(pos, cam, las, high_las, mod_cyc, ncycs_ignore) #note: the warmup cycles are named differently but sequence is same as other fast pos cycles
# 
#	add_sequence!(high_las, "laser high", trues(nsamps_cycle)) #will append one of these for every cycle as we go
#	#Now iterate through slices, lags, and trials
#	for s = 1:length(slice_zs)
#		fwd_ctr, back_ctr = find_bidi_locations(mean_cyc, samprate(pos), slice_zs[s])
#		for l = 1:Base.length(lags)
#			fwd_exp_name = "fwd_exp_$l"
#			fwd_flash_name = "fwd_flash_$l"
#			back_exp_name = "back_exp_$l"
#			back_flash_name = "back_flash_$l"
#			#A lag of 0.0s would mean that the laser pulse is centered on the measured stack location
#			#Positive lags mean that the pulse centers precede the measured stack location
#			fwd_exp, fwd_flash = get_trial(fwd_ctr, lags[l], samprate(pos), nsamps_exp, nsamps_flash, nsamps_sweep, true)
#			back_exp, back_flash = get_trial(back_ctr, lags[l], samprate(pos), nsamps_exp, nsamps_flash, nsamps_sweep, false)
#			add_sequence!(cam, back_exp_name, back_exp)
#			add_sequence!(las, back_flash_name, back_flash)
#			add_sequence!(cam, fwd_exp_name, fwd_exp)
#			add_sequence!(las, fwd_flash_name, fwd_flash)
#			for t = 1:ntrials
#				append!(pos, pos_cyc_name)
#				append!(cam, fwd_exp_name)
#				append!(cam, back_exp_name)
#				append!(las, fwd_flash_name)
#				append!(las, back_flash_name)
#				append!(high_las, "laser_high")
#			end
#		end
#	end
#    if nsamps_offset > 0 #add another positioner cycle because we have too many cam and las samples
#        append!(pos, pos_cyc_name)
#        append!(cam, "extra_cyc", fill(false, nsamps_cycle-nsamps_offset))
#        append!(las, "extra_cyc")
#        append!(high_las, "extra_cychigh", fill(true, nsamps_cycle-nsamps_offset))
#    end
#    return [pos; cam; all_las...]
#end
#
#function slicetimings_step2_gen(pos::T1, mon_cyc, las_name::AbstractString, cam_name::AbstractString, lags, slice_zs) where {T1<:ImagineSignal}
#    mod_cyc = get_samples(pos, first(sequence_names(pos))) #relies on the cycle being the first (and only) subsequence
#    @assert all(sequence_names(pos).==pos_cyc_name)
#    rig_sigs = rigtemplate(rig_name(pos); sample_rate = samprate(pos))
#    cam = getname(rig_sigs, cam_name)
#    las = high_las = all_las = -1
#    @assert rig == "ocpi-2"
#    if las_name != "all lasers"
#        las = getname(rig_sigs, "all lasers")
#        high_las = getname(rig_sigs, las_name)
#        all_las = [las;high_las]
#    else
#        error("Please specify a specific laser, not all lasers")
#    end
#    if !ImagineInterface.ispos(pos) error("Only positioner signals are currently supported.") end
#    if !ImagineInterface.iscam(cam) error("$cam_name is not a camera trigger signal.") end
#    if !ImagineInterface.islas(las) error("$las_name is not a laser trigger signal.") end
#
#    slicetimings_step2_gen!(deepcopy(pos), mod_cyc, mon_cyc, las, las_high, cam, lags, slice_zs) where {T1<:ImagineSignal}
#end
#
#
