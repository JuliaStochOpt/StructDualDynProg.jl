# Similar to JuMP/test/solvers.jl

function try_import(name::Symbol)
    try
        @eval import $name
        return true
    catch e
        return false
    end
end

#using ECOS
#solver = ECOS.ECOSSolver(verbose=false)
clp = try_import(:Clp)
isclp(solver) = contains(string(typeof(solver)),"ClpSolver")
glp = try_import(:GLPKMathProgInterface)
gur = try_import(:Gurobi)

lp_solvers = Any[]
clp && push!(lp_solvers, Clp.ClpSolver())
glp && push!(lp_solvers, GLPKMathProgInterface.GLPKSolverLP())
gur && push!(lp_solvers, Gurobi.GurobiSolver())
