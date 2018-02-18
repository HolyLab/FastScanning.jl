#slicetiming_experiment runs a pilot piezo-only recording, reads piezo monitor samples, and returns an imagine procedure with the following attributes:
#The SignalGenerator generates command signals with two phases:
#   Phase 1:  A slow (0.2Hz) stack is acquired with camera and laser pulses timed to match target z locations determined by z_spacing.
#   Phase 2:  Many piezo cycles with two camera and laser pulses per cycle.  The piezo cycle is taken from the recs input.
#           Each pulse is timed at a different lag relative to the target z location determined by z_spacing (and used to place slices in Phase 1)
#The process function expects the recorded signals from the above commands as inputs, and it outputs time indices indicating when to time slices so that z location
#when imaging fast matches the location of images from the template stack of phase 1 above.
#Optionally pass a Calibration to be applied to the piezo MON signal before doing all this (if, for example, you know the constant lag in the MON output)

default_toffsets() = [0.0s:0.00005s:0.0022s...]

function slicetiming_experiment(out_bname::AbstractString, rig::AbstractString, pos_name, las_name, cam_name, pstart, pstop, stack_rate, z_spacing, z_pad, ncycs_mean=ceil(Int, 20.0/ustrip(inv(stack_rate))); sample_rate = 100000Hz, cal=-1, toffsets = default_toffsets())
    mod_cyc = vcat(gen_bidi_pos(pstart, pstop, 1/stack_rate, sample_rate)...)
    return slicetiming_experiment(out_bname, rig, pos_name, las_name, cam_name, mod_cyc, z_spacing, z_pad; ncycs_mean=ncycs_mean, sample_rate=sample_rate, cal=cal, toffsets=toffsets)
end

function slicetiming_experiment(out_bname::AbstractString, rig::AbstractString, pos_name, las_name, cam_name, mod_cyc, z_spacing, z_pad; ncycs_mean=ceil(Int, 20.0/ustrip(length(mod_cyc)*sample_rate)), sample_rate = 100000Hz, cal=-1, toffsets = default_toffsets())
    coms0 = pos_commands(rig, pos_name, mod_cyc, ncycs_mean; sample_rate =sample_rate)
    pos = getpositioners(coms0)
    pos_mon = getpositionermonitors(coms0)
    #run commands
    warn("This method must be run on the microscope computer while the piezo is on and connected (both MON and MOD connections)")
    write_commands(out_bname * "_piezo_only.json", [pos; pos_mon], 0, 0, 0.01s; exp_trig_mode = ["External Start"], isbidi=true, skip_validation=true)
    recs = Imagine.run_imagine(out_bname * "_piezo_only", [pos; pos_mon]; ai_trig_dest = "PFI2", ao_trig_dest = "PFI1", trigger_source = "Port2/Line0", skip_validation = true)
    return slicetiming_experiment([pos; pos_mon], recs, las_name, cam_name, z_spacing, z_pad; cal = cal, toffsets=toffsets)
end

function slicetiming_experiment(coms, recs, las_name::AbstractString, cam_name::AbstractString, z_spacing::HasLengthUnits, z_pad::HasLengthUnits; cal = -1, toffsets = default_toffsets())
    pos, pos_mon = pos_mod_mon(coms, recs)
    if cal == -1
        pos_mon_cal = pos_mon
    else
        pos_mon_cal = apply(cal, pos_mon)
    end
    mon_samps = get_samples(pos_mon_cal)
    nsamps_cycle = full_length(first(sequences(pos)))
    ncycs_ignore = ceil(Int, (ustrip(samprate(pos))*5) / nsamps_cycle) #5 seconds of warmup
    mean_cyc = mean_cycle(mon_samps, nsamps_cycle; start_idx = ncycs_ignore*nsamps_cycle+1)
    #set pmin and pmax empirically
    pmin = minimum(mean_cyc)
    pmax = maximum(mean_cyc)
    slice_zs = ImagineInterface.slice_positions(pmin, pmax, z_spacing, z_pad)
    slicetiming_experiment(pos, mean_cyc, las_name, cam_name, slice_zs; toffsets=toffsets)
end

function slicetiming_experiment(pos::ImagineSignal, mon_cyc, las_name::AbstractString, cam_name::AbstractString, slice_zs; toffsets = default_toffsets())
    gen_desc = "Generate a set camera and laser commands to find slice timings empirically given a completed dynamic positioner recording"
    sig_gen = SignalGenerator(gen_desc, slicetiming_commands, (pos, mon_cyc, las_name, cam_name, toffsets, slice_zs))
    desc = "Returns slice timings to align forward and reverse stacks at the desired spacing for high-speed volume acquisitions"
    analyze_f = (coms, img) -> (slicetimings(img, toffsets, mon_cyc, length(slice_zs))..., slice_zs) #returns forward toffsets, backward toffsets, and the z locations of slices
    #analyze_f = (coms, img) -> calc_lag(coms, img, toffsets_pre, ncycs_ignore, mon_cyc, nsamps_offset)
    return ImagineProcedure(desc, sig_gen, analyze_f)
end

#Similar to slicetiming_experiment but only checks one slice per stack, and uses the result to estimate the lag between the piezo's true position and the MON signal
#function pos_mon_lag_experiment(pos::ImagineSignal, mon_cyc::ImagineSignal, las_name::AbstractString, cam_name::AbstractString, slice_zs)

function pos_mon_lag_experiment(pos::ImagineSignal, mon_cyc, las_name::AbstractString, cam_name::AbstractString)
    #set pmin and pmax empirically
    pmin = minimum(mon_cyc)
    pmax = maximum(mon_cyc)
    pctr = [(pmax-pmin)/2 * pmin]
    @show toffsets = default_toffsets() #These worked to find the lag on the digital piezo of OCPI2
    gen_desc = "Generate a set camera and laser commands to find piezo monitor lag empirically given a completed dynamic positioner recording"
    sig_gen = SignalGenerator(gen_desc, slicetiming_commands, (pos, mon_cyc, las_name, cam_name, toffsets, [pctr]))
    desc = "Returns mean of two measurements: the monitor lag when taking an image while sweeping forward and also while sweeping back"
    analyze_f = (coms, img) -> mean(slicetimings(img, toffsets, mon_cyc, 1)) #take the mean of the forward and reverse lag (hopefully they are about equal)
    return ImagineProcedure(desc, sig_gen, analyze_f)
end

pos_mon_lag_experiment(coms, recs, las_name::AbstractString, cam_name::AbstractString) = pos_mon_lag_experiment(pos_mod_mon(coms,recs)..., las_name, cam_name)

function pos_mon_lag_experiment(pos::ImagineSignal, pos_mon::ImagineSignal, las_name::AbstractString, cam_name::AbstractString)
    mon_samps = get_samples(pos_mon)
    nsamps_cycle = full_length(first(sequences(pos)))
    ncycs_ignore = ceil(Int, (ustrip(samprate(pos))*5) / nsamps_cycle) #5 seconds of warmup
    mon_cyc = mean_cycle(mon_samps, nsamps_cycle; start_idx = ncycs_ignore*nsamps_cycle+1)
    pos_mon_lag_experiment(pos, mon_cyc, las_name, cam_name)
end
