module SDDP

using TimerOutputs

using StructDualDynProg
const SOI = StructDualDynProg.StochOptInterface

# SDDP algorithm for StochOptInterface

# Utils
include("comp.jl")

# Path Sampler
include("sampler.jl")

include("path.jl")
include("sddp.jl")

end
