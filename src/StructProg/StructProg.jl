module StructProg

using Compat, Compat.LinearAlgebra, Compat.SparseArrays

using TimerOutputs

using StochOptInterface
const SOI = StochOptInterface
using StructDualDynProg

using JuMP
using StructJuMP

using CutPruners

# Utils
include("comp.jl")

# Cut Generator
include("cutgen.jl")

# Generic implementation of Stochastic Program
include("cutstore.jl")
include("solver.jl")
include("nlds.jl")
include("graph.jl")

# Cut
include("cut.jl")

# StructJuMP interface
include("interface.jl")

end
