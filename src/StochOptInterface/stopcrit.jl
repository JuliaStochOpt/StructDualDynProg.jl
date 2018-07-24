abstract type AbstractStoppingCriterion end

"""
    stop(s::AbstractStoppingCriterion, info)

Determine whether the SDDP algorithm should stop using the information stored in `info`.
"""
function stop(s::AbstractStoppingCriterion, info::Info)
    error("`stop' function not defined for $(typeof(s))")
end

"""
$(TYPEDEF)

Stops if `lhs` *or* `rhs` want to stop.
"""
mutable struct OrStoppingCriterion <: AbstractStoppingCriterion
    lhs::AbstractStoppingCriterion
    rhs::AbstractStoppingCriterion
end

function stop(s::OrStoppingCriterion, info::Info)
    stop(s.lhs, info) || stop(s.rhs, info)
end

function Base.:(|)(lhs::AbstractStoppingCriterion, rhs::AbstractStoppingCriterion)
    OrStoppingCriterion(lhs, rhs)
end

"""
$(TYPEDEF)

Stops if `lhs` *and* `rhs` want to stop.
"""
mutable struct AndStoppingCriterion <: AbstractStoppingCriterion
    lhs::AbstractStoppingCriterion
    rhs::AbstractStoppingCriterion
end

function stop(s::AndStoppingCriterion, info::Info)
    stop(s.lhs, info) && stop(s.rhs, info)
end

function Base.:(&)(lhs::AbstractStoppingCriterion, rhs::AbstractStoppingCriterion)
    AndStoppingCriterion(lhs, rhs)
end

"""
$(TYPEDEF)

Stops if `iter` ≧ `limit`.
"""
mutable struct IterLimit <: AbstractStoppingCriterion
    limit::Int
end

function stop(s::IterLimit, info::Info)
    niterations(info) >= s.limit
end

"""
$(TYPEDEF)

Stops if there was less than or equal to `limit` cuts added in the iteration.
For instance, `CutLimit(0)` stops when there are no cuts added.
"""
mutable struct CutLimit <: AbstractStoppingCriterion
    limit::Int
end

function stop(s::CutLimit, info::Info)
    niterations(info) > 0 && nfcuts(last_timer(info)) + nocuts(last_timer(info)) <= s.limit
end


"""
$(TYPEDEF)

Stops if total time of execution is greater than the time limit specified.
For instance, `TimeLimit(100)` stops after 100s.
"""
mutable struct TimeLimit <: AbstractStoppingCriterion
    timelimit::Float64
end

function stop(s::TimeLimit, info::Info)
    niterations(info) > 0 && total_time(info) > s.timelimit
end


"""
$(TYPEDEF)

Stops if `z_UB - α * σ/√K - tol < z_LB < z_UB + α * σ/√K + tol` and `σ / √K > β * max(1, |z_LB|))`
"""
mutable struct Pereira <: AbstractStoppingCriterion
    α::Float64
    β::Float64
    tol::Float64

    Pereira(α=2.0, β=0.05, tol=1e-6) = new(α, β, tol)
end

function stop(s::Pereira, info::Info)
    if niterations(info) > 0
        result = last_result(info)
        z_UB = result.upperbound
        z_LB = result.lowerbound
        K = result.npaths
        σ = result.σ_UB

        @assert K >= 0
        σ1 = σ / √K
        # On the test optimize_stock with Clp, z_LB = -2, z_UB = -1.999999999999 and σ1 = 0
        # this shows the necessicity for a tolerance
        σ2 = s.α * σ1 + s.tol
        z_UB - σ2 <= z_LB <= z_UB + σ2 && σ1 < s.β * max(1, abs(z_LB))
    else
        false
    end
end
