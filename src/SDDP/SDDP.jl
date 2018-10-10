module SDDP

using Compat, Compat.LinearAlgebra

using TimerOutputs

# SDDP algorithm for StochOptInterface
using StochOptInterface
const SOI = StochOptInterface

using StructDualDynProg

# Utils
include("comp.jl")

# Path Sampler
include("sampler.jl")

include("path.jl")
include("algorithm.jl")

end
