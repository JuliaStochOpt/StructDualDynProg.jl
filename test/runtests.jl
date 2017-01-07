using JuMP
using StructJuMP
using StochasticDualDynamicProgramming
using Base.Test

# solver independent tests
include("cutmanager.jl")

# load a solver
include("solvers.jl")

# solver dependent tests
include("optimize_stock.jl")
include("prob5.2.jl")
