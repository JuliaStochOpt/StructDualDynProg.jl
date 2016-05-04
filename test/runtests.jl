#using GLPKMathProgInterface
#using ECOS
using Clp
using JuMP
using StructJuMP
using StochasticDualDynamicProgramming
using Base.Test

#solver = ECOS.ECOSSolver(verbose=false)
solver = Clp.ClpSolver()
#solver = GLPKSolverLP()

include("prob5.2.jl")
include("prob5.2_3stages.jl")
