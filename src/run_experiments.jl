#This code only gets loaded when the functionality of Imagine.jl is available (on a Windows machine with NIDAQmx installed)

using Imagine

function slicetiming_experiment(out_bname::AbstractString, rig::AbstractString, pos_name, las_name, cam_name, mod_cyc, z_spacing, z_pad; ncycs_mean=ceil(Int, 20.0/ustrip(length(mod_cyc)*sample_rate)), sample_rate = 100000Hz, cal=-1, toffsets = default_toffsets(), allow_shifts=true, allow_rotations=false, subpixel=true)
    coms0 = pos_commands(rig, pos_name, mod_cyc, ncycs_mean; sample_rate =sample_rate)
    pos = getpositioners(coms0)
    pos_mon = getpositionermonitors(coms0)
    #run commands
    warn("This method must be run on the microscope computer while the piezo is on and connected (both MON and MOD connections)")
    write_commands(out_bname * "_piezo_only.json", [pos; pos_mon], 0, 0, 0.01s; exp_trig_mode = ["External Start"], isbidi=true, skip_validation=true)
    recs = Imagine.run_imagine(out_bname * "_piezo_only", [pos; pos_mon]; ai_trig_dest = "PFI2", ao_trig_dest = "PFI1", trigger_source = "Port2/Line0", skip_validation = true)
    return slicetiming_experiment([pos; pos_mon], recs, las_name, cam_name, z_spacing, z_pad; cal = cal, toffsets=toffsets, allow_shifts=allow_shifts, allow_rotations=allow_rotations, subpixel=subpixel)
end
