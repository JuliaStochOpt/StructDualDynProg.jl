__precompile__()

module StructDualDynProg

using DocStringExtensions

# TODO This module should be replaced by the package https://github.com/JuliaStochOpt/StochOptInterface.jl
include("StochOptInterface/StochOptInterface.jl")
const SOI = StochOptInterface

# Generic implementation of Stochastic Program supporting the conversion from a StructJuMP model to a Stochastic Program
include("StructProg/StructProg.jl")
export StructProg

# SDDP algorithm on top of these abstract components
include("SDDP/SDDP.jl")
export SDDP

# Wait and see value computation for StructProg
include("WaitAndSee/waitandsee.jl")
export WaitAndSee

end # module
