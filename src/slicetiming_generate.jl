#Generate an ouput ImagineSignal for commanding the piezo to repeat the sample sequence `one_cyc` `ncycles` times.
#Appends those samples to `pos` and also returns dummy camera and laser signals so that Imagine doesn't complain
#`ncycles` default to 20 about seconds-worth of cycles based on the sample rate of `pos`
function pos_commands(pos::ImagineSignal, one_cyc, ncycles=ceil(Int, 20.0 / ustrip(inv(samprate(pos))*length(one_cyc))))
    @assert isempty(pos)
    rig_sigs = rigtemplate(rig_name(pos); sample_rate = samprate(pos))
    add_sequence!(pos, "bidi_cycle", one_cyc)
    cycle_nsamps = length(one_cyc)
    cam = first(getcameras(rig_sigs)) #not used
    las = first(getlasers(rig_sigs)) #not used
    add_sequence!(cam, "quiet", falses(cycle_nsamps)) #just to pass validation
    for i = 1:ncycles
        append!(pos, "bidi_cycle")
        append!(cam, "quiet")
        append!(las, "quiet")
    end
    pos_mon = getname(rig_sigs, monitor_name(pos))
    return [pos;cam;las;pos_mon]
end

function pos_commands(rig::AbstractString, pos_name::AbstractString, pstart, pstop, stack_rate, ncycles=ceil(Int, 20.0/ustrip(inv(stack_rate))); sample_rate = 100000Hz)
    mod_cyc = vcat(gen_bidi_pos(pstart, pstop, 1/stack_rate, sample_rate)...)
    return pos_commands(rig, pos_name, mod_cyc, ncycles; sample_rate = sample_rate)
end

function pos_commands(rig::AbstractString, pos_name::AbstractString, mod_cyc, ncycles=ceil(Int, 20.0/ustrip(inv(stack_rate))); sample_rate = 100000Hz)
    rig_sigs = rigtemplate(rig; sample_rate = sample_rate)
    pos = getname(rig_sigs, pos_name)
    return pos_commands(pos, mod_cyc, ncycles)
end

function append_warmup!(pos, cam, las, high_las, pos_cycle, ncycs_warmup)
    nsamps_cyc = Base.length(pos_cycle) #should always be even
    add_sequence!(cam, "cam_warmup", falses(nsamps_cyc))
    add_sequence!(las, "las_warmup", falses(nsamps_cyc))
    add_sequence!(pos, "pos_warmup", pos_cycle)
    if rig_name(pos) == "ocpi-2"
        add_sequence!(high_las, "high_las_warmup", trues(nsamps_cyc))
    end
    for i = 1:ncycs_warmup
        append!(pos, "pos_warmup")
        append!(cam, "cam_warmup")
        append!(las, "las_warmup")
        append!(high_las, "high_las_warmup")
    end
	return pos, cam, las, high_las
end

function get_trial_intervals(ctr_i, lag, sample_rate, nsamps_exp, nsamps_flash)
	nsamps_from_ctr = div(nsamps_flash,2)
	nsamps_lag = ImagineInterface.calc_num_samps(lag, sample_rate)
    flash_start = ctr_i-nsamps_lag-nsamps_from_ctr
    flash_stop = ctr_i-nsamps_lag+nsamps_from_ctr
    exp_stop = flash_stop
    exp_start = flash_stop-nsamps_exp+1
    flash_itv = ClosedInterval(flash_start, flash_stop)
    exp_itv = ClosedInterval(exp_start, exp_stop)
    return exp_itv, flash_itv
end

function check_overlap(itvs::Vector{T}) where {T<:ClosedInterval}
    for i = 2:length(itvs)
        if !isempty(intersect(itvs[i-1], itvs[i]))
            error("Pulse intervals overlap")
        end
    end
    return itvs
end

#given a measured piezo cycle, a set of desired slice locations `slice_zs`, and a set of forward and reverse lags,
#return a tuple (cam_samps, las_samps, shift) with sample vectors for controlling camera exposure and laser pulsing
function get_cyc_pulses(mon_cyc, slice_zs, fwd_lags, back_lags, exp_time, flash_time, sample_rate::HasInverseTimeUnits)
	nsamps_exp = ImagineInterface.calc_num_samps(exp_time, sample_rate)
	nsamps_flash = ImagineInterface.calc_num_samps(flash_time, sample_rate)
    @assert iseven(length(mon_cyc))
    @assert issorted(slice_zs) || issorted(slice_zs; rev=true)
    nsamps_sweep = div(length(mon_cyc),2)
    cams_fwd = ClosedInterval{Int}[]
    cams_back = ClosedInterval{Int}[]
    lasers_fwd = ClosedInterval{Int}[]
    lasers_back = ClosedInterval{Int}[]
    nslices = length(slice_zs)
    #get_trial_intervals can return intervals lying outside the cycle bounds (negative first indices and last indices greater than 2*nsamps_sweep)
    #check extrema of intervals (smin and smax below)
    #if both of these extrema lie outside of 2*nsamps_sweep, throw an error
    #if only one of these does, it's okay.  return cam_samps, las_samps, and the shift value
    for s = 1:nslices
        ctr_fwd, ctr_back = find_bidi_locations(mon_cyc, sample_rate, slice_zs[s])
        cam_fwd, las_fwd = get_trial_intervals(ctr_fwd, fwd_lags[s], sample_rate, nsamps_exp, nsamps_flash)
        cam_back, las_back = get_trial_intervals(ctr_back, back_lags[s], sample_rate, nsamps_exp, nsamps_flash)
        #cam_back+=nsamps_sweep
        #las_back+=nsamps_sweep
        push!(cams_fwd, cam_fwd)
        push!(cams_back, cam_back)
        push!(lasers_fwd, las_fwd)
        push!(lasers_back, las_back)
    end
    las_itvs = vcat(lasers_fwd, lasers_back)
    cam_itvs = vcat(cams_fwd, cams_back)
    shift = 0
    all_itvs = vcat(las_itvs, cam_itvs)
    smin = minimum(minimum.(all_itvs))
    smax = maximum(maximum.(all_itvs))
    if smin < 1 && smax > 2*nsamps_sweep
        error("Cannot create pulse sequence without overlapping pulses across cycles")
    elseif !(smin >0 && (smax <= 2*nsamps_sweep))
        if abs(smin-1) > abs(smax - 2*nsamps_sweep)
            @show shift = smin-1
        else
            @show shift = smax - 2*nsamps_sweep
        end
    end
    las_samps = ImagineInterface.gen_pulses(2*nsamps_sweep, check_overlap(las_itvs.-shift))
    cam_samps = ImagineInterface.gen_pulses(2*nsamps_sweep, check_overlap(cam_itvs.-shift))
    return cam_samps, las_samps, shift #shift is the shift that should be applied to sample vectors to get them into proper alignment
end

function append_template!(pos, cam, las, high_las, slice_zs, exp_time, flash_time; freq=0.2Hz) #appends a stack sequence to each signal input. Images only on the forward sweep.
    #pad by 5μm if we have room
    z_min = minimum(slice_zs)
    z_max = maximum(slice_zs)
    z_pad = minimum([abs(z_min - 0.0μm), abs(800μm - z_max), 5.0μm])
    z_min-=z_pad
    z_max+=z_pad
    mod_cyc = vcat(gen_bidi_pos(z_min, z_max, 1/freq, samprate(pos))...)
    z_range = z_max - z_min #should also match the command (not necessarily the monitored signal)
    pos_speed = freq * z_range #in microns per second
    pad_nsecs = z_pad / pos_speed #nsecs to cover the z_pad distance
    pad_nsamps = ImagineInterface.calc_num_samps(pad_nsecs, samprate(pos))
    flash_cyc, cam_cyc, monshft = flash_cam_cycs(mod_cyc, mod_cyc, slice_zs, pad_nsamps, flash_time, exp_time, samprate(pos))
    #only take pictures on the forward sweep
    midpoint = div(length(flash_cyc),2) + 1
    flash_cyc[midpoint:end] = false
    cam_cyc[midpoint:end] = false
    append!(pos, "pos_template", mod_cyc)
    append!(cam, "cam_template", cam_cyc)
    append!(las, "las_template", flash_cyc)
    append!(high_las, "highlas_template", trues(length(mod_cyc)))
end

function multistack_bidi(cam_name, las_name, pos_name, pos_cyc, mon_cyc, nstacks_fwd, slice_zs, fwd_toffsets, back_toffsets, exp_time, flash_time, sample_rate, ncycs_warmup)
	ocpi2 = rigtemplate("ocpi-2"; sample_rate = sample_rate)
	pos_new = getname(ocpi2, pos_name)
	cam_new = getname(ocpi2, cam_name)
	las_new = getname(ocpi2, "all lasers")
	high_las_new = getname(ocpi2, las_name)

	EmpiricalTiming.append_warmup!(pos_new, cam_new, las_new, high_las_new, pos_cyc, ncycs_warmup)

    @show nsamps_offset = ImagineAnalyses.mon_delay(pos_cyc, mon_cyc)
	#temporal shift to make calculations easier (will shift things back later)
    mon_cyc = circshift(mon_cyc, -nsamps_offset)
	if nsamps_offset < 0
		error("can't handle negative offset")
	elseif nsamps_offset > 0
		append!(las_new, "delay", fill(false, nsamps_offset))
		append!(cam_new, "delay")
		append!(high_las_new, "delayhigh", fill(true, nsamps_offset))
	end
	cam_samps, las_samps, shift = get_cyc_pulses(mon_cyc, slice_zs, fwd_toffsets, back_toffsets, exp_time, flash_time, sample_rate)
	@assert shift == 0 #TODO: handle this

	nsamps_cyc = length(pos_cyc)
	pos_cyc_name = "pos_cyc"
	add_sequence!(pos_new, pos_cyc_name, pos_cyc)
	add_sequence!(cam_new, "cam_cyc", cam_samps)
	add_sequence!(las_new, "las_cyc", las_samps)
	add_sequence!(high_las_new, "high_cyc", trues(nsamps_cyc))

	for i = 1:nstacks_fwd
		append!(pos_new, "pos_cyc")
		append!(cam_new, "cam_cyc")
		append!(las_new, "las_cyc")
		append!(high_las_new, "high_cyc")
	end
    if nsamps_offset > 0 #add another positioner cycle because we have too many cam and las samples
        append!(pos_new, pos_cyc_name)
        append!(cam_new, "extra_cyc", fill(false, nsamps_cyc-nsamps_offset))
        append!(las_new, "extra_cyc")
        append!(high_las_new, "extra_cychigh", fill(true, nsamps_cyc-nsamps_offset))
    end
	return [pos_new; cam_new; las_new; high_las_new]
end

function slicetiming_commands!(pos::T1, mod_cyc, mon_cyc, las, high_las, cam, lags, slice_zs) where {T1<:ImagineSignal}
    nsamps_cyc = Base.length(mod_cyc) #should always be even
    nsamps_sweep = div(nsamps_cyc,2)
    ncycs_ignore = ceil(Int, (ustrip(samprate(pos))*5) / nsamps_cyc) #5 seconds of warmup
    #@show nsamps_offset = mon_lag_nsamps(mod_cyc, mon_cyc)
    @show nsamps_offset = ImagineAnalyses.mon_delay(mod_cyc, mon_cyc)
	#temporal shift to make calculations easier (will shift things back later)
    mon_cyc = circshift(mon_cyc, -nsamps_offset)
    ntrials=1 #number of trials per lag value per z slice
    #We want to deliver the laser pulse during the global shutter time (tglobal in PCO camera manual)
    #If the user images with full chip then the global exposure segment doesn't begin until _10ms_ after the exposure starts.
    #Therefore with a 1ms laser pulse we need an 11ms exposure time, and the pulse is placed in the last 1ms 
    exp_time = 0.0020s
    flash_time = 0.0005s
    empty!(pos)
	if nsamps_offset < 0
		error("Cannot handle negative offset (but could with a bit more coding work)")
	elseif nsamps_offset > 0
		append!(las, "delay", fill(false, nsamps_offset))
		append!(cam, "delay")
		append!(high_las, "delayhigh", fill(true, nsamps_offset))
	end
    nsamps_flash = ImagineInterface.calc_num_samps(flash_time, samprate(pos))
    if iseven(nsamps_flash)
	    nsamps_flash -= 1 #force odd so that fwd and reverse flashes are easier to align
    end
	nsamps_exp = ImagineInterface.calc_num_samps(exp_time, samprate(pos))

	#Sequence from here:
	# one slow image stack
	# several fast "warmup" piezo cycles with no camera / laser activity (ncycs_ignore)
	# ntrials * nlags piezo cycles in which we take a pair of slices for each cycle
	# repeat above step for each slice, so a total of ntrials * nlags * nslices cycles (and twice that many images)
    ncycs_ignore = ceil(Int, (ustrip(samprate(pos))*5) / nsamps_cyc) #5 seconds of warmup
	append_template!(pos, cam, las, high_las, slice_zs, exp_time, flash_time) #appends slow stack sequence to each signal input
	append_warmup!(pos, cam, las, high_las, mod_cyc, ncycs_ignore) #note: the warmup cycles are named differently but sequence is same as other fast pos cycles
    
    pos_cyc_name = "fast_pos_cyc"
	add_sequence!(pos, pos_cyc_name, mod_cyc)
	add_sequence!(high_las, "laser_high", trues(nsamps_cyc)) #will append one of these for every cycle as we go
	#Now iterate through slices, lags, and trials
	for s = 1:length(slice_zs)
		fwd_ctr, back_ctr = find_bidi_locations(mon_cyc, samprate(pos), slice_zs[s])
		for l = 1:Base.length(lags)
            cyc_exp_name = "cyc_exp_$(l)_slice_$s"
            cyc_flash_name = "cyc_flash_$(l)_slice_$s"
            has_partial = false #do we need to add a partial cycle before and after each cycle due to exposure pulse timing?
            cam_samps, las_samps, timeshift = get_cyc_pulses(mon_cyc, [slice_zs[s]], [lags[l];], [lags[l];], exp_time, flash_time, samprate(pos))
			add_sequence!(cam, cyc_exp_name, cam_samps)
			add_sequence!(las, cyc_flash_name, las_samps)
            #if timeshift is nonzero then the imaging pulses lie partially outside the cycle and we need an extra cycle
            if timeshift < 0
                warn("Adding an extra piezo cycle in slice $s with lag $l for exposures or flashes that begin before the cycle")
                has_partial = true
                add_sequence!(cam, "cam_partial_before_$(s)_$(l)", falses(nsamps_cyc+timeshift))
                add_sequence!(las, "las_partial_before_$(s)_$(l)", falses(nsamps_cyc+timeshift))
                add_sequence!(cam, "cam_partial_after_$(s)_$(l)",  falses(-timeshift))
                add_sequence!(las, "las_partial_after_$(s)_$(l)",  falses(-timeshift))
                @show s
                @show l
            elseif timeshift > 0
                has_partial = true
                warn("Adding an extra piezo cycle in slice $s with lag $l for exposures or flashes that end after the cycle")
                add_sequence!(cam, "cam_partial_before_$(s)_$(l)", falses(timeshift))
                add_sequence!(las, "las_partial_before_$(s)_$(l)", falses(timeshift))
                add_sequence!(cam, "cam_partial_after_$(s)_$(l)", falses(nsamps_cyc-timeshift))
                add_sequence!(las, "las_partial_after_$(s)_$(l)", falses(nsamps_cyc-timeshift))
                @show s
                @show l
            end
            if has_partial
                append!(pos, pos_cyc_name)
                append!(cam, "cam_partial_before_$(s)_$(l)")
                append!(las, "las_partial_before_$(s)_$(l)")
                append!(high_las, "laser_high")
            end
			for t = 1:ntrials
				append!(pos, pos_cyc_name)
				append!(cam, cyc_exp_name)
				append!(las, cyc_flash_name)
				append!(high_las, "laser_high")
			end
            if has_partial
                append!(cam, "cam_partial_after_$(s)_$(l)")
                append!(las, "las_partial_after_$(s)_$(l)")
            end
		end
	end
    if nsamps_offset > 0 #add another positioner cycle because we have too many cam and las samples
        append!(pos, pos_cyc_name)
        append!(cam, "extra_cyc", fill(false, nsamps_cyc-nsamps_offset))
        append!(las, "extra_cyc")
        append!(high_las, "extra_cychigh", fill(true, nsamps_cyc-nsamps_offset))
    end
    @assert length(pos) == length(cam) == length(las) == length(high_las)
    return [pos; cam; las; high_las]
end

function slicetiming_commands(pos::T1, mon_cyc, las_name::AbstractString, cam_name::AbstractString, lags, slice_zs) where {T1<:ImagineSignal}
    mod_cyc = get_samples(pos, first(sequence_names(pos))) #relies on the cycle being the first (and only) subsequence
    @assert all(sequence_names(pos).==first(sequence_names(pos)))
    rig_sigs = rigtemplate(rig_name(pos); sample_rate = samprate(pos))
    cam = getname(rig_sigs, cam_name)
    las = high_las = all_las = -1
    @assert rig_name(pos) == "ocpi-2"
    if las_name != "all lasers"
        las = getname(rig_sigs, "all lasers")
        high_las = getname(rig_sigs, las_name)
        all_las = [las;high_las]
    else
        error("Please specify a specific laser, not all lasers")
    end
    if !ImagineInterface.ispos(pos) error("Only positioner signals are currently supported.") end
    if !ImagineInterface.iscam(cam) error("$cam_name is not a camera trigger signal.") end
    if !ImagineInterface.islas(las) error("$las_name is not a laser trigger signal.") end

    slicetiming_commands!(deepcopy(pos), mod_cyc, mon_cyc, las, high_las, cam, lags, slice_zs)
end
