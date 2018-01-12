module EmpiricalTiming

using Unitful, ImagineInterface, ImagineAnalyses, ImagineProcedures, Images, AffineTransforms, ProgressMeter
using BlockRegistrationScheduler #for tweaking imperfect 2D alignment of slices

import Unitful:s, Hz, μm, V

include("slicetiming_analysis.jl")
include("bsearch.jl")
include("slicetiming_generate.jl")
include("experiments.jl")

export get_cyc_pulses, pos_mon_lag_experiment, slicetiming_experiment

end # module
