#Generate an ouput ImagineSignal for commanding the piezo to repeat the sample sequence `one_cyc` `ncycles` times.
#Appends those samples to `pos` and also returns dummy camera and laser signals so that Imagine doesn't complain
#`ncycles` default to 20 about seconds-worth of cycles based on the sample rate of `pos`
function pos_commands(pos::ImagineSignal, one_cyc::Vector, ncycles=ceil(Int, 20.0 / ustrip(inv(samprate(pos))*length(one_cyc))))
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
    return [pos;cam;las]
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
    nsamps_cycle = Base.length(pos_cycle) #should always be even
    add_sequence!(cam, "cam_warmup", falses(nsamps_cycle))
    add_sequence!(las, "las_warmup", falses(nsamps_cycle))
    add_sequence!(pos, "pos_warmup", pos_cycle)
    if rig_name(pos) == "ocpi-2"
        add_sequence!(high_las, "high_las_warmup", trues(nsamps_cycle))
    end
    for i = 1:ncycs_warmup
        append!(pos, "pos_warmup")
        append!(cam, "cam_warmup")
        append!(las, "las_warmup")
        append!(high_las, "high_las_warmup")
    end
	return pos, cam, las, high_las
end

function get_trial_intervals(ctr_i, lag, sample_rate, nsamps_exp, nsamps_flash, nsamps_sweep, is_fwd::Bool)
	nsamps_from_ctr = div(nsamps_flash,2)
	nsamps_lag = ImagineInterface.calc_num_samps(lag, sample_rate)
#    @show sample_rate
#    @show nsamps_lag
#    @show lag
#    @show ctr_i
#    @show nsamps_exp
#    @show nsamps_flash
    shift = 0
	if is_fwd
        flash_start = ctr_i-nsamps_lag-nsamps_from_ctr
        flash_stop = ctr_i-nsamps_lag+nsamps_from_ctr
        exp_start = ctr_i-nsamps_lag+nsamps_from_ctr-nsamps_exp-1
        exp_stop = ctr_i-nsamps_lag+nsamps_from_ctr
        if exp_start < 1 #exposure always starts before flash so only need to check that.
            shift = exp_start - 1
            flash_start -= shift
            flash_stop -= shift
            exp_start -= shift
            exp_stop -= shift
        end
	else
        flash_start = ctr_i-nsamps_lag-nsamps_from_ctr-nsamps_sweep
        flash_stop = ctr_i-nsamps_lag+nsamps_from_ctr-nsamps_sweep
        exp_start = ctr_i-nsamps_lag+nsamps_from_ctr-nsamps_exp-1-nsamps_sweep
        exp_stop = ctr_i-nsamps_lag+nsamps_from_ctr-nsamps_sweep
        if exp_stop > 2*nsamps_sweep #exposure always stops at the same time as flash so only need to check that.
            shift = 2*nsamps_sweep - exp_stop
            flash_start -= shift
            flash_stop -= shift
            exp_start -= shift
            exp_stop -= shift
        end
	end
    flash_itv = ClosedInterval(flash_start, flash_stop)
    exp_itv = ClosedInterval(exp_start, exp_stop)
    return exp_itv, flash_itv, shift
end

function get_trial(ctr_i, lag, sample_rate, nsamps_exp, nsamps_flash, nsamps_sweep, is_fwd::Bool)
    exp_itv, flash_itv, shift = get_trial_intervals(ctr_i, lag, sample_rate, nsamps_exp, nsamps_flash, nsamps_sweep, is_fwd)
	las_samps = ImagineInterface.gen_pulses(nsamps_sweep, [flash_itv;])
	cam_samps = ImagineInterface.gen_pulses(nsamps_sweep, [exp_itv;])
	return cam_samps, las_samps, shift
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
    append!(high_las, "highlas_template", falses(length(mod_cyc)))
end

function slicetiming_commands!(pos::T1, mod_cyc, mon_cyc, las, high_las, cam, lags, slice_zs) where {T1<:ImagineSignal}
    nsamps_cycle = Base.length(mod_cyc) #should always be even
    nsamps_sweep = div(nsamps_cycle,2)
    ncycs_ignore = ceil(Int, (ustrip(samprate(pos))*5) / nsamps_cycle) #5 seconds of warmup
    #@show nsamps_offset = mon_lag_nsamps(mod_cyc, mon_cyc)
    @show nsamps_offset = ImagineAnalyses.mon_delay(mod_cyc, mon_cyc)
	#temporal shift to make calculations easier (will shift things back later)
    mon_cyc = circshift(mon_cyc, -nsamps_offset)
    ntrials=1 #number of trials per lag value per z slice
    #We want to deliver the laser pulse during the global shutter time (tglobal in PCO camera manual)
    #If the user images with full chip then the global exposure segment doesn't begin until _10ms_ after the exposure starts.
    #Therefore with a 1ms laser pulse we need an 11ms exposure time, and the pulse is placed in the last 1ms 
    exp_time = 0.0030s
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
    ncycs_ignore = ceil(Int, (ustrip(samprate(pos))*5) / nsamps_cycle) #5 seconds of warmup
	append_template!(pos, cam, las, high_las, slice_zs, exp_time, flash_time) #appends slow stack sequence to each signal input
	append_warmup!(pos, cam, las, high_las, mod_cyc, ncycs_ignore) #note: the warmup cycles are named differently but sequence is same as other fast pos cycles
    
    pos_cyc_name = "fast_pos_cyc"
	add_sequence!(pos, pos_cyc_name, mod_cyc)
	add_sequence!(high_las, "laser_high", trues(nsamps_cycle)) #will append one of these for every cycle as we go
	#Now iterate through slices, lags, and trials
	for s = 1:length(slice_zs)
		fwd_ctr, back_ctr = find_bidi_locations(mon_cyc, samprate(pos), slice_zs[s])
		for l = 1:Base.length(lags)
            fwd_exp_name = "fwd_exp_$(l)_slice_$s"
            fwd_flash_name = "fwd_flash_$(l)_slice_$s"
            back_exp_name = "back_exp_$(l)_slice_$s"
            back_flash_name = "back_flash_$(l)_slice_$s"
			#A lag of 0.0s would mean that the laser pulse is centered on the measured stack location
			#Positive lags mean that the pulse centers precede the measured stack location
			fwd_exp, fwd_flash, timeshift_fwd = get_trial(fwd_ctr, lags[l], samprate(pos), nsamps_exp, nsamps_flash, nsamps_sweep, true)
			back_exp, back_flash, timeshift_back = get_trial(back_ctr, lags[l], samprate(pos), nsamps_exp, nsamps_flash, nsamps_sweep, false)
			add_sequence!(cam, back_exp_name, back_exp)
			add_sequence!(las, back_flash_name, back_flash)
			add_sequence!(cam, fwd_exp_name, fwd_exp)
			add_sequence!(las, fwd_flash_name, fwd_flash)
            #if timeshifts are nonzero then the imaging pulses lie partially outside the cycle and we need an extra cycle or two.
            if timeshift_fwd != 0
                error("timeshift_fwd = $timeshift_fwd")
            elseif timeshift_back != 0
                error("timeshift_back = $timeshift_back")
            end
#            if timeshift_fwd < 0
#                append!(pos, pos_cyc_name)
#                append!(cam, "overflow_fwd_lag_$l", cam[)
#            end
			for t = 1:ntrials
				append!(pos, pos_cyc_name)
				append!(cam, fwd_exp_name)
				append!(cam, back_exp_name)
				append!(las, fwd_flash_name)
				append!(las, back_flash_name)
				append!(high_las, "laser_high")
			end
		end
	end
    if nsamps_offset > 0 #add another positioner cycle because we have too many cam and las samples
        append!(pos, pos_cyc_name)
        append!(cam, "extra_cyc", fill(false, nsamps_cycle-nsamps_offset))
        append!(las, "extra_cyc")
        append!(high_las, "extra_cychigh", fill(true, nsamps_cycle-nsamps_offset))
    end
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
