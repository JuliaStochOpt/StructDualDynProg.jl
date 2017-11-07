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
