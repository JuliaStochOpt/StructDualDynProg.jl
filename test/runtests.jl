using ECOS
using JuMP
using StructJuMP
using StochasticDualDynamicProgramming
using Base.Test

solver = ECOS.ECOSSolver(verbose=false)

include("prob5.2.jl")
include("prob5.2_3stages.jl")
