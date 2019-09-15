#using Unitful, ImagineInterface, ImagineAnalyses, Images, ProgressMeter
#import Unitful:s, μm, V
#
##calculate lag mapping the piezo mon signal to the hardware response
##Procedure:
##   The idea is to use the camera as an independent measurement device instead of the piezo sensor signal
##   Delivere several paired-exposure pulse trials with bead sample:
##        first pulse always occurs in the middle of the forward stroke of the piezo (bidirectional motion)
##        second pulse always occurs in the middle of the reverse stroke of the piezo, offset by a set of lag times.
##        ssd on pairs of images to determine at what lag time ssd is minimized.  That is the monitor lag.
#
##lags is which lag values to check
#function mon_lag_sigs0(rig::AbstractString, pos_name::AbstractString, pos_start = 540μm, pos_stop = 630μm, stack_rate = 10*inv(Unitful.s))
#    sample_rate = 100000 * inv(Unitful.s)
#    rig_sigs = rigtemplate(rig; sample_rate = sample_rate)
#    posfwd, posback = gen_bidi_pos(pos_start, pos_stop, 1/stack_rate, sample_rate)
#    pos = getname(rig_sigs, pos_name)
#    add_sequence!(pos, "bidi_cycle", vcat(posfwd, posback))
#    cycle_nsamps = length(posfwd) + length(posback)
#    cam = first(getcameras(rig_sigs))
#    las = first(getlasers(rig_sigs))
#    add_sequence!(cam, "quiet", falses(cycle_nsamps)) #just to pass validation
#    ncycles = 300
#    for i = 1:ncycles
#        append!(pos, "bidi_cycle")
#        append!(cam, "quiet")
#        append!(las, "quiet")
#    end
#    return [pos;cam;las]
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
#function gen_mon_lag_procedure(coms, recs, las_name::AbstractString, cam_name::AbstractString)
#    pos = getpositioners(coms)
#    if length(pos) > 1
#        error("Provide only one positioner command.")
#    end
#    pos = pos[1]
#    pos_mon = getpositionermonitors(recs)
#    if length(pos_mon) > 1
#        error("Provide only one positioner monitor.")
#    end
#    pos_mon = pos_mon[1]
#    nsamps_cycle = full_length(first(sequences(pos)))
#    ncycs_ignore = 50
#    ntrials=10
#    @show nsamps_offset = mon_lag_nsamps(pos, pos_mon, ncycs_ignore)
#    mon_samps = get_samples(pos_mon)
#    mean_cyc = mean_cycle(mon_samps, nsamps_cycle; start_idx = ncycs_ignore*nsamps_cycle+1)
#    mean_cyc = circshift(mean_cyc, -nsamps_offset)
#    #@show lags = [0.0002s:0.00005s:0.0007s...]  #These worked to find the lag on the analog piezo of OCPI2
#    @show lags_pre = [0.0005s:0.00005s:0.002s...] #These worked to find the lag on the digital piezo of OCPI2
#    gen_desc = "Generate a set camera and laser commands to measure positioner monitor lag given a completed dynamic positioner recording"
#    sig_gen = SignalGenerator(gen_desc, mon_lag_sigs, (pos, mean_cyc, las_name, cam_name, lags_pre, ncycs_ignore, ntrials, nsamps_offset))
#    rig = rig_name(pos)
#    mon_name = name(pos_mon)
#    desc = "Calculates lag mapping the $rig $mon_name sensor signal to its true hardware response.  Returns a Calibration."
#    analyze_f = (coms, img) -> calc_lag(coms, img, lags_pre, ncycs_ignore, mean_cyc, nsamps_offset)
#    return ImagineProcedure(desc, sig_gen, analyze_f)
#end
#
#function calc_lag(coms, img, lags, ncycs_ignore::Int, mean_cyc, nsamps_offset::Int)
#    #mean_cyc has already been shifted by offset so that it aligns well with the command
#    sampr = samprate(first(coms))
#    fwd_ctr, back_ctr = find_bidi_centers(mean_cyc, sampr)
#    lag_pre = select_lag_pre(coms, img, lags, ncycs_ignore)
#    #NOTE: the lag chosen here isn't the signal lag, it needs more processing to find that.
#    #The lag chosen here is the temporal displacement from the MON measurement at the exposure moment to the moment when the
#    #forward and reverse sweeps cross the same location (that crossing takes place somewhere between the corrected z location
#    #and the original target z location
#    #so the true value is roughly half of this measurement
#    #a more precise way to find the true value is to advance two pointers until they reach the same z value:
#    #the first pointer begins at the fwd flash location and steps forward in time, the second
#    #pointer begins at the corrected second flash location and steps forward in time.
#    #The (equal) time interval traveled by the pointers when they meet is the lag
#    fwd_i = fwd_ctr
#    back_i = back_ctr - ImagineInterface.calc_num_samps(lag_pre, sampr)
#    lag_nsamps = 0
#    while back_i < length(mean_cyc)
#        if mean_cyc[fwd_i] < mean_cyc[back_i]
#            lag_nsamps+=1
#            fwd_i+=1
#            back_i+=1
#        else
#            return Calibration(upreferred(lag_nsamps / sampr), 0.0V, 1.0)
#        end
#    end
#    error("Couldn't find a lag value that makes sense.  Check inputs.")
#end
#
#
##Analyzes pairs of images to find the lag (in lags) at which image difference is minimized
##The lags vector is assumed to be ordered in the same way that trials were ordered when acquiring imgs
#function select_lag_pre(coms, img, lags, ncycs_ignore)
#    pos = first(getpositioners(coms))
#    ncycs = length(sequence_names(pos))
#    npairs = size(img,3) / 2
#    if !isinteger(npairs)
#        error("Expected an odd number of images (each trial requires an image pair).")
#    elseif npairs < 1
#        error("Not enough cycles.")
#    else
#        npairs = Int(npairs)
#    end
#    if abs(ncycs_ignore - abs(ncycs - npairs)) > 1 #allow a difference of one because sometimes an extra cycle gets added when the MON signal lags a lot
#        error("The number of cycles in the command signal disagrees with the image count and ncycs_ignore")
#    end
#    nlags = length(lags)
#    ntrials = npairs/nlags
#    if !isinteger(ntrials)
#        error("The number of image pairs should be evenly divisible by the number of lag conditions.")
#    else
#        ntrials = Int(ntrials)
#    end
#    sz_spatial = size(img)[1:2]
#    img = reshape(img, sz_spatial..., 2, ntrials, nlags)
#    ssds = zeros(ntrials, nlags)
#    #TODO: should also do a rigid 2D registration when comparing images because there is lateral wobble in fwd vs reverse images
#    @showprogress for l = 1:nlags
#        for t = 1:ntrials
#            ssds[t,l] = ssd(img[:,:,1,t,l], img[:,:,2,t,l])
#        end
#    end
#    @show mean_ssds = mean(ssds, 1)[:]
#    return lags[indmin(mean_ssds)]
#end
#
#function find_bidi_centers(cyc::AbstractVector, sampr::HasInverseTimeUnits)
#    @show pmin = minimum(cyc)
#    @show pmax = maximum(cyc)
#    prange = pmax - pmin
#    pctr = pmin + prange/2 #we will place the flash in the center
#    stack_rate = 2 * (sampr / length(cyc))
#    pos_avg_speed = stack_rate * prange #in microns per second
#    pad_secs = (1/stack_rate)/2 #wait half a cycle before accepting threshold crossings (helps with noise)
#    pad_nsamps = round(Int, pad_secs * sampr)
#    mon_ctr_is = ImagineAnalyses.find_circular(cyc, [pctr;], pad_nsamps)[1]
#    @assert Base.length(mon_ctr_is) == 2 #should have just two crossings, one fwd one back
#    fwd_ctr = mon_ctr_is[1]
#    back_ctr = mon_ctr_is[2]
#    return fwd_ctr, back_ctr
#end
#
##ncycs_ignore determines how many to skip before starting exposure pulses
##...and also how many cycles of pos_mon to skip when calculating the mean cycle waveform
#function mon_lag_sigs(pos::T1, mean_cyc, las_name::AbstractString, cam_name::AbstractString, lags, ncycs_ignore, ntrials, nsamps_offset::Int) where {T1<:ImagineSignal}
#    exp_time = 0.011s
#    flash_time = 0.001s
#    pos = deepcopy(pos)
#    rig = rig_name(pos)
#    #We want to deliver the laser pulse during the global shutter time (tglobal in PCO camera manual)
#    #If the user images with full chip then the global exposure segment doesn't begin until _10ms_ after the exposure starts.
#    #Therefore with a 1ms laser pulse we need an 11ms exposure time, and the pulse is placed in the last 1ms
#    rig_sigs = rigtemplate(rig; sample_rate = samprate(pos))
#    cam = getname(rig_sigs, cam_name)
#    las = -1
#    high_las = -1
#    all_las = -1
#    if rig == "ocpi-2"
#        if las_name != "all lasers"
#            las = getname(rig_sigs, "all lasers")
#            high_las = getname(rig_sigs, las_name)
#            all_las = [las;high_las]
#        else
#            error("Please specify a specific laser, not all lasers")
#        end
#    else
#        las = getname(rig_sigs, las_name)
#        all_las = [las;]
#    end
#
#    if !ImagineInterface.ispos(pos)
#        error("Only positioner signals are currently supported.")
#    end
#    if !ImagineInterface.iscam(cam)
#        error("$cam_name is not a camera trigger signal.")
#    end
#    if !ImagineInterface.islas(las)
#        error("$las_name is not a laser trigger signal.")
#    end
#    
#    pos_cyc_name = first(sequence_names(pos))
#    pos_cycle = get_samples(pos, pos_cyc_name)
#    @assert all(sequence_names(pos).==pos_cyc_name)
#    empty!(pos)
#
#    nsamps_cycle = Base.length(pos_cycle) #should always be even
#    add_sequence!(cam, "wait_to_start", falses(nsamps_cycle))
#    if rig == "ocpi-2"
#        add_sequence!(high_las, "laser_high", trues(nsamps_cycle))
#    end
#    for i = 1:ncycs_ignore
#        append!(pos, pos_cyc_name)
#        append!(cam, "wait_to_start")
#        append!(las, "wait_to_start")
#        append!(high_las, "laser_high")
#    end
#    #After the ignored cycles each trial has this structure:
#    #pos: pos_cyc_name
#    #cam: "fwd_exp" -> "back_exp_$i" where i is an index into the lags vector
#    #las: "fwd_flash" -> "back_flash_$i"
#    #in order to know when/where to place the camera and laser pulses we must look at the pos_mon signal
#
#    fwd_ctr, back_ctr = find_bidi_centers(mean_cyc, samprate(pos))
#    
#    nsamps_flash = ImagineInterface.calc_num_samps(flash_time, samprate(pos))
#    if iseven(nsamps_flash)
#	nsamps_flash -= 1 #force odd so that fwd and reverse flashes are easier to align
#    end
#    nsamps_from_ctr = div(nsamps_flash,2)
#    #forward flash and exposure are always timed the same, so we can calculate them just once
#    flash_itv_fwd = ClosedInterval(fwd_ctr-nsamps_from_ctr, fwd_ctr+nsamps_from_ctr) 
#    nsamps_exp = ImagineInterface.calc_num_samps(exp_time, samprate(pos))
#    #This aligns the end of the exposure with the end of the flash (allows global exposure)
#    exp_itv_fwd = ClosedInterval(fwd_ctr+nsamps_from_ctr-nsamps_exp-1, fwd_ctr+nsamps_from_ctr)
#    nsamps_sweep = div(nsamps_cycle,2)
#    fwd_flash = ImagineInterface.gen_pulses(nsamps_sweep, [flash_itv_fwd;])
#    fwd_exp = ImagineInterface.gen_pulses(nsamps_sweep, [exp_itv_fwd;])
#    if nsamps_offset < 0
#        error("Cannot handle negative monitor lag")
#    elseif nsamps_offset > 0
#        append!(las, "delay", fill(false, nsamps_offset))
#        append!(cam, "delay")
#        append!(high_las, "delayhigh", fill(true, nsamps_offset))
#    end
#    add_sequence!(las, "fwd_flash", fwd_flash)
#    add_sequence!(cam, "fwd_exp", fwd_exp)
#    for l = 1:Base.length(lags)
#        #calculate flash and exposure samples for backward sweep
#        back_exp_name = "back_exp_$l"
#        back_flash_name = "back_flash_$l"
#        #A lag of 0.0s would mean that the laser pulse is centered on the measured stack center (back_ctr)
#        #Positive lags mean that the pulse centers precede back_ctr
#        nsamps_lag = ImagineInterface.calc_num_samps(lags[l], samprate(pos))
#        flash_itv_back = ClosedInterval(back_ctr-nsamps_lag-nsamps_from_ctr-nsamps_sweep, back_ctr-nsamps_lag+nsamps_from_ctr-nsamps_sweep)
#        exp_itv_back = ClosedInterval(back_ctr-nsamps_lag+nsamps_from_ctr-nsamps_exp-1-nsamps_sweep, back_ctr-nsamps_lag+nsamps_from_ctr-nsamps_sweep)
#        back_flash = ImagineInterface.gen_pulses(nsamps_sweep, [flash_itv_back;])
#        back_exp = ImagineInterface.gen_pulses(nsamps_sweep, [exp_itv_back;])
#        add_sequence!(cam, back_exp_name, back_exp)
#        add_sequence!(las, back_flash_name, back_flash)
#        for t = 1:ntrials
#            append!(pos, pos_cyc_name)
#            append!(cam, "fwd_exp")
#            append!(cam, back_exp_name)
#            append!(las, "fwd_flash")
#            append!(las, back_flash_name)
#            if rig == "ocpi-2"
#                append!(high_las, "laser_high")
#            end
#        end
#    end
#    if nsamps_offset > 0 #add another positioner cycle because we have too many cam and las samples
#        append!(pos, pos_cyc_name)
#        append!(cam, "extra_cyc", fill(false, nsamps_cycle-nsamps_offset))
#        append!(las, "extra_cyc")
#        append!(high_las, "extra_cychigh", fill(true, nsamps_cycle-nsamps_offset))
#    end
#    return [pos; cam; all_las...]
#end
#
##create procedures for each positioner on each rig
#rigs = ImagineInterface.RIGS
#rigs = filter(x->x!="dummy-6002", rigs)
#for rig in rigs
#    rig_sigs = rigtemplate(rig)
#    pos_mons = getpositionermonitors(rig_sigs)
#    for mon in pos_mons
#        mon_name = name(mon)
#        mod_name = ImagineInterface.actuator_name(mon)
#        desc = "Produces a procedure for measuring a lag mapping the $rig $mon_name sensor signal to its true hardware response.  Returns an ImagineProcedure."
#        sig_gen = SignalGenerator("Generate piezo command on which monitor lag calculation will be based.", mon_lag_sigs0, (rig, mod_name))
#        analyze_f = gen_mon_lag_procedure
#        p = ImagineProcedure(desc, sig_gen, analyze_f)
#        #now store it
#        store_calibration_procedure(rig, daq_channel(mon), p)
#    end
#end
