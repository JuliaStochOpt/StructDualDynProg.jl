# Solution at one state, different than Solution with is the full solution
abstract type AbstractSolution end

"""
    feasibility_cut(sol::AbstractSolution)

Returns the tuple `(a, α)` representing the feasibility cut ``⟨a, x⟩ ≧ α`` certified by this solution.
"""
function feasibility_cut end

"""
    optimality_cut(sol::AbstractSolution)

Returns the tuple `(a, α)` representing the optimality cut ``⟨a, x⟩ + θ ≧ α`` certified by this solution.
"""
function optimality_cut end

"""
    getstatus(sol::AbstractSolution)

Returns the status of the solution `sol`.
"""
function getstatus end

"""
    getobjectivevalue(sol::AbstractSolution)

Returns the objective value of the solution `sol` *including* the part of the objective depending on `θ`.
"""
function getobjectivevalue end

"""
    getstateobjectivevalue(sol::AbstractSolution)

Returns the objective value of the solution `sol` *excluding* the part of the objective depending on `θ`.
"""
function getstateobjectivevalue end


"""
    getstatevalue(sol::AbstractSolution)

Returns the value of the state of solution `sol`.
"""
function getstatevalue end

"""
    getθvalue(sp::AbstractStochasticProgram, tr::AbstractTransition, sol::AbstractSolution)

Returns the value of the θ in the solution `sol` of node `SOI.get(sp, SOI.Source(), tr)` for its transition `tr`.
This assumes that `node` is using `MultiCutGenerator`.

    getθvalue(sp::AbstractStochasticProgram, node, sol::AbstractSolution)

Returns the value of the θ in the solution `sol` of node `node`.
This assumes that `node` is using `AvgCutGenerator`.
"""
function getθvalue end

abstract type AbstractSolutionPool end

"""
    allfeasible(pool::AbstractSolutionPool)

Return a `Bool` indicating whether all transitions current solved were feasible.
"""
function allfeasible end

"""
    allbounded(pool::AbstractSolutionPool)

Return a `Bool` indicating whether all transitions current solved were bounded.
"""
function allbounded end

"""
    hassolution(pool::AbstractSolutionPool, tr::AbstractTransition)

Return a `Bool` indicating whether the solution pool `pool` has a solution for transition `tr`.
"""
function hassolution end

"""
    getsolution(pool::AbstractSolutionPool)

Return the solution for the source of all the transitions.

    getsolution(pool::AbstractSolutionPool, tr::AbstractTransition)

Return the solution for transition `tr` in the solution pool `pool`.
"""
function getsolution end
