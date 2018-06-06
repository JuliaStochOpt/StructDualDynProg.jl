module SDDP

using StructDualDynProg
const SOI = StructDualDynProg.StochOptInterface

# SDDP algorithm for StochOptInterface

# Path Sampler
include("sampler.jl")

include("path.jl")
include("sddp.jl")

end
