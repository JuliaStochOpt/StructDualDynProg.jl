# Similar to JuMP/test/solvers.jl

function try_import(name::Symbol)
    try
        @eval import $name
        return true
    catch e
        return false
    end
end

grb = false && try_import(:Gurobi)
isgrb(solver) = occursin("GurobiSolver", string(typeof(solver)))
cpx = false && try_import(:CPLEX)
iscpx(solver) = occursin("CplexSolver", string(typeof(solver)))
xpr = false && try_import(:Xpress)
clp = false && try_import(:Clp)
isclp(solver) = occursin("ClpSolver", string(typeof(solver)))
glp = try_import(:GLPK)
msk = false && try_import(:Mosek)
ismsk(solver) = occursin("MosekSolver", string(typeof(solver)))

lp_solvers = Any[]
grb && push!(lp_solvers, optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0))
cpx && push!(lp_solvers, optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_SCRIND" => 0, "CPX_PARAM_REDUCE" => 0))
xpr && push!(lp_solvers, optimizer_with_attributes(Xpress.Optimizer, "OUTPUTLOG" => 0))
clp && push!(lp_solvers, Clp.Optimizer)
glp && push!(lp_solvers, GLPK.Optimizer)
msk && push!(lp_solvers, optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true))
