module StructProg

using DocStringExtensions

using StructDualDynProg
const SOI = StructDualDynProg.StochOptInterface

using JuMP
using StructJuMP

using CutPruners

# Generic implementation of Stochastic Program
include("cutstore.jl")
include("solver.jl")
include("nlds.jl")
include("graph.jl")

# StructJuMP interface
include("interface.jl")

end
