using FastScanning
using Test
using Unitful, AxisArrays, Images, ImagineInterface, ImagineProcedures
import Unitful:s, Hz, μm

ocpi2 = rigtemplate("ocpi-2")
pos_mon = getname(ocpi2, "axial piezo monitor")
pstart = 15.0μm
pstop = 115.0μm
stack_rate = 10.0Hz
z_spacing = 5.0μm
z_pad = 5.0μm
coms = FastScanning.pos_commands("ocpi-2", "axial piezo", pstart, pstop, stack_rate; sample_rate = 100000Hz)
pos = getpositioners(coms)[1]
fake_samps = get_samples(pos)
append!(pos_mon, "fake_samps", fake_samps) #perfect with no lag
pr = slicetiming_experiment(coms, [pos_mon;], "488nm laser", "camera1", z_spacing, z_pad; subpixel=false)
st_coms = outputs(pr)

lags = [0.0000s:0.0001s:0.002s...] #note: must match what's in experiments.jl
ntrials = 1

nconditions = length(lags)
slice_zs = ImagineInterface.slice_positions(pstart, pstop, z_spacing, z_pad)
nslices = length(slice_zs)
best_conds_fwd = fill(1, nslices)
best_conds_back = fill(2, nslices)
imgs = fake_slicetiming_run(nconditions, ntrials, nslices, best_conds_fwd, best_conds_back)

fwd_lags, back_lags, slice_zs2 = process(pr, imgs)
#@assert slice_zs2.==slice_zs #not equal due to finite sample rate
@test lags[best_conds_fwd] == fwd_lags
@test lags[best_conds_back] == back_lags


#### again with subpixel == true
pr = slicetiming_experiment(coms, [pos_mon;], "488nm laser", "camera1", z_spacing, z_pad; subpixel=true)
st_coms = outputs(pr)
imgs = fake_slicetiming_run(nconditions, ntrials, nslices, best_conds_fwd, best_conds_back)

fwd_lags, back_lags, slice_zs2 = process(pr, imgs)
@test lags[best_conds_fwd] == fwd_lags
@test lags[best_conds_back] == back_lags
