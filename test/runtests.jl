using JuMP
using CutPruners
using StructJuMP
using StochasticDualDynamicProgramming
using Base.Test

# solver independent tests
include("stopcrit.jl")
include("paths.jl")

# load a solver
include("solvers.jl")

# solver dependent tests
include("optimize_stock.jl")
include("prob5.2.jl")
