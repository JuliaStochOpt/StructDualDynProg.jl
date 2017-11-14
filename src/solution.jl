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
    getθvalue(sp::AbstractStochasticProgram, node, child, sol::AbstractSolution)

Returns the value of the θ in the solution `sol` of node `node` for its child `child`.
This assumes that `node` is using `MultiCutGenerator`.

    getθvalue(sp::AbstractStochasticProgram, node, sol::AbstractSolution)

Returns the value of the θ in the solution `sol` of node `node`.
This assumes that `node` is using `AvgCutGenerator`.
"""
function getθvalue end
