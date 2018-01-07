module EmpiricalTiming

using Unitful, ImagineInterface, ImagineAnalyses, ImagineProcedures, Images, ProgressMeter
using BlockRegistrationScheduler #for tweaking imperfect 2D alignment of slices

import Unitful:s, Hz, μm, V

include("slicetiming_analysis.jl")
include("slicetiming_generate.jl")
include("experiments.jl")

export pos_mon_lag_experiment, slicetiming_experiment

end # module
