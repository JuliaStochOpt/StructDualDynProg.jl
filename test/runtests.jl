using JuMP
using StructJuMP
using StochasticDualDynamicProgramming
using Base.Test

#using ECOS
#solver = ECOS.ECOSSolver(verbose=false)
#using Clp
#solver = Clp.ClpSolver()
using Gurobi
solver = Gurobi.GurobiSolver(OutputFlag=0)
#using GLPKMathProgInterface
#solver = GLPKSolverLP()

include("prob5.2.jl")
include("prob5.2_3stages.jl")
include("prob5.2_3stages_serial.jl")
