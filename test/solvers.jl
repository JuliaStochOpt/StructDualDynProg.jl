# Similar to JuMP/test/solvers.jl

function try_import(name::Symbol)
    try
        @eval import $name
        return true
    catch e
        return false
    end
end

grb = try_import(:Gurobi)
isgrb(solver) = contains(string(typeof(solver)),"GurobiSolver")
cpx = try_import(:CPLEX)
iscpx(solver) = contains(string(typeof(solver)),"CplexSolver")
xpr = try_import(:Xpress)
clp = try_import(:Clp)
isclp(solver) = contains(string(typeof(solver)),"ClpSolver")
glp = try_import(:GLPKMathProgInterface)
msk = false && try_import(:Mosek)
ismsk(solver) = contains(string(typeof(solver)),"MosekSolver")

lp_solvers = Any[]
grb && push!(lp_solvers, Gurobi.GurobiSolver(OutputFlag=0))
cpx && push!(lp_solvers, CPLEX.CplexSolver(CPX_PARAM_SCRIND=0, CPX_PARAM_REDUCE=0))
xpr && push!(lp_solvers, Xpress.XpressSolver(OUTPUTLOG=0))
clp && push!(lp_solvers, Clp.ClpSolver())
glp && push!(lp_solvers, GLPKMathProgInterface.GLPKSolverLP())
msk && push!(lp_solvers, Mosek.MosekSolver(QUIET=true))
