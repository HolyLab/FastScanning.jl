#reads piezo monitor samples from recs and returns an imagine procedure with the following attributes:
#The SignalGenerator generates command signals with two phases:
#   Phase 1:  A slow (0.2Hz) stack is acquired with camera and laser pulses timed to match target z locations determined by z_spacing.
#   Phase 2:  Many piezo cycles with two camera and laser pulses per cycle.  The piezo cycle is taken from the recs input.
#           Each pulse is timed at a different lag relative to the target z location determined by z_spacing (and used to place slices in Phase 1)
#The process function expects the recorded signals from the above commands as inputs, and it outputs time indices indicating when to time slices so that z location
#when imaging fast matches the location of images from the template stack of phase 1 above.
#Optionally pass a Calibration to be applied to the piezo MON signal before doing all this (if, for example, you know the constant lag in the MON output)
#TODO:make another method that runs the pilot positioner command and then passes the recorded MON signal in the `recs` argument
#   The signature for that method:
#   function slicetiming_experiment(pos_name, las_name::AbstractString, cam_name::AbstractString, freq, z_start, z_stop, z_spacing::HasLengthUnits, z_pad::HasLengthUnits; cal = Calibration())
function slicetiming_experiment(coms, recs, las_name::AbstractString, cam_name::AbstractString, z_spacing::HasLengthUnits, z_pad::HasLengthUnits; cal = -1)
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
    slicetiming_experiment(pos, mean_cyc, las_name, cam_name, slice_zs)
end

function slicetiming_experiment(pos::ImagineSignal, mon_cyc, las_name::AbstractString, cam_name::AbstractString, slice_zs)
    #@show lags = [0.0002s:0.00005s:0.0007s...]  #These worked to find the lag on the analog piezo of OCPI2
    @show lags = [0.0005s:0.00005s:0.0012s...] #These worked to find the lag on the digital piezo of OCPI2
    gen_desc = "Generate a set camera and laser commands to find slice timings empirically given a completed dynamic positioner recording"
    sig_gen = SignalGenerator(gen_desc, slicetiming_commands, (pos, mon_cyc, las_name, cam_name, lags, slice_zs))
    desc = "Returns slice timings to align forward and reverse stacks at the desired spacing for high-speed volume acquisitions"
    analyze_f = (comss, img) -> (slicetimings(img, lags, mon_cyc, length(slice_zs))..., slice_zs) #returns forward lags, backward lags, and the z locations of slices
    #analyze_f = (coms, img) -> calc_lag(coms, img, lags_pre, ncycs_ignore, mon_cyc, nsamps_offset)
    return ImagineProcedure(desc, sig_gen, analyze_f)
end

#Similar to slicetiming_experiment but only checks one slice per stack, and uses the result to estimate the lag between the piezo's true position and the MON signal
#function pos_mon_lag_experiment(pos::ImagineSignal, mon_cyc::ImagineSignal, las_name::AbstractString, cam_name::AbstractString, slice_zs)

function pos_mon_lag_experiment(pos::ImagineSignal, mon_cyc, las_name::AbstractString, cam_name::AbstractString)
    #set pmin and pmax empirically
    pmin = minimum(mon_cyc)
    pmax = maximum(mon_cyc)
    pctr = [(pmax-pmin)/2 * pmin]
    #@show lags = [0.0002s:0.00005s:0.0007s...]  #These worked to find the lag on the analog piezo of OCPI2
    @show lags = [0.0005s:0.00005s:0.002s...] #These worked to find the lag on the digital piezo of OCPI2
    gen_desc = "Generate a set camera and laser commands to find piezo monitor lag empirically given a completed dynamic positioner recording"
    sig_gen = SignalGenerator(gen_desc, slicetiming_commands, (pos, mon_cyc, las_name, cam_name, lags, [pctr]))
    desc = "Returns mean of two measurements: the monitor lag when taking an image while sweeping forward and also while sweeping back"
    analyze_f = (coms, img) -> mean(slicetimings(img, lags, mon_cyc, 1)) #take the mean of the forward and reverse lag (hopefully they are about equal)
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
