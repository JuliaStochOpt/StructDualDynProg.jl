# π, σ and ρ do not really make sense alone so only
# their product will T, h, d, e is stored
mutable struct Solution
    status::Symbol
    objval
    objvalx
    objvalxuray
    x
    xuray # unbouded ray
    θ
    θuray # unbouded ray
    πT
    πh
    σd
    ρe
    function Solution(status::Symbol, objval)
        new(status, objval, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing)
    end
end

# Feasibility cut
# D = π T
# d = π h + σ d
function getfeasibilitycut(sol::Solution)
    (sol.πT, sol.πh + sol.σd)
end

# Optimality cut
# E = π T
# e = π h + ρ e + σ d
function getoptimalitycut(sol::Solution)
    (sol.πT, sol.πh + sol.σd + sol.ρe)
end
