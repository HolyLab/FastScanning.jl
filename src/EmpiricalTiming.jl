module EmpiricalTiming

using Unitful, ImagineInterface, ImagineAnalyses, ImagineProcedures
using Images, CoordinateTransformations, ProgressMeter, FixedPointNumbers, IntervalSets
using Statistics, Distributed
using RegisterQD, RegisterMismatch
using Requires

import Unitful:s, Hz, μm, V

include("slicetiming_analysis.jl")
#include("bsearch.jl")
include("slicetiming_generate.jl")
include("experiments.jl")

function __init__()
    @require Imagine="28fce3f0-1dbc-11e9-10e3-e71e6ceb4e7c" include("run_experiments.jl")
end

export default_toffsets, pos_commands, get_cyc_pulses, pos_mon_lag_experiment, slicetiming_experiment

end # module
