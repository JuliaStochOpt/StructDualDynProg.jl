module StructProg

using DocStringExtensions

using StructDualDynProg
const SOI = StructDualDynProg.StochOptInterface

using TimerOutputs

using JuMP
using StructJuMP

using CutPruners

# Utils
include("comp.jl")

# Cut
include("cut.jl")
# Cut Generator
include("cutgen.jl")

# Generic implementation of Stochastic Program
include("cutstore.jl")
include("solver.jl")
include("nlds.jl")
include("graph.jl")

# StructJuMP interface
include("interface.jl")

end
