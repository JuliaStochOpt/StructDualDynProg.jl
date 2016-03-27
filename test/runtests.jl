using Gurobi
using ECOS
using JuMP
using StructJuMP
using StochasticDualDynamicProgramming
using Base.Test

misocp_solver = GurobiSolver(OutputFlag=0)
socp_solver = ECOS.ECOSSolver(verbose=false)

include("prob5.2.jl")
#include("prob5.2_3stages.jl")
