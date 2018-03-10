module EmpiricalTiming

using Unitful, Imagine, ImagineInterface, ImagineAnalyses, ImagineProcedures, Images, AffineTransforms, CoordinateTransformations, ProgressMeter, FixedPointNumbers
using BlockRegistrationScheduler #for tweaking imperfect 2D alignment of slices
using BlockRegistration, RegisterOptimize, RegisterMismatch

import Unitful:s, Hz, μm, V

include("slicetiming_analysis.jl")
include("bsearch.jl")
include("slicetiming_generate.jl")
include("experiments.jl")

export default_toffsets, get_cyc_pulses, pos_mon_lag_experiment, slicetiming_experiment

end # module
