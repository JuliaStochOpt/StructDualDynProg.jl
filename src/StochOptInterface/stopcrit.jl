abstract type AbstractStoppingCriterion end

"""
    stop(s::AbstractStoppingCriterion, to::TimerOutput, stats::AbstractStats, totalstats::AbstractStats)

Returns whether the SDDP algorithm should stop.
If `totalstats.niterations` is 0, no iteration has already been done, otherwise, the `niterations`th iteration has just finished.
This iteration used `stats.npaths` paths and generated `nfcuts(to)` (resp. `nocuts(to)`) new feasibility (resp. optimality) cuts.
The lower bound is now `totalstats.lowerbound` and the upper bound has mean `totalstats.upperbound` and variance `totalstats.σ_UB`.
"""
function stop(s::AbstractStoppingCriterion, to::TimerOutput, stats::AbstractStats, totalstats::AbstractStats)
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

function stop(s::OrStoppingCriterion, to::TimerOutput, stats::AbstractStats, totalstats::AbstractStats)
    stop(s.lhs, to, stats, totalstats) || stop(s.rhs, to, stats, totalstats)
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

function stop(s::AndStoppingCriterion, to::TimerOutput, stats::AbstractStats, totalstats::AbstractStats)
    stop(s.lhs, to, stats, totalstats) && stop(s.rhs, to, stats, totalstats)
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

function stop(s::IterLimit, to::TimerOutput, stats::AbstractStats, totalstats::AbstractStats)
    totalstats.niterations >= s.limit
end

"""
$(TYPEDEF)

Stops if there was less than or equal to `limit` cuts added in the iteration.
For instance, `CutLimit(0)` stops when there are no cuts added.
"""
mutable struct CutLimit <: AbstractStoppingCriterion
    limit::Int
end

function stop(s::CutLimit, to::TimerOutput, stats::AbstractStats, totalstats::AbstractStats)
    totalstats.niterations > 0 && nfcuts(to) + nocuts(to) <= s.limit
end


"""
$(TYPEDEF)

Stops if total time of execution is greater than the time limit specified.
For instance, `TimeLimit(100)` stops after 100s.
"""
mutable struct TimeLimit <: AbstractStoppingCriterion
    timelimit::Float64
end

function stop(s::TimeLimit, to::TimerOutput, stats::AbstractStats, totalstats::AbstractStats)
    totalstats.niterations > 0 && totalstats.time > s.timelimit
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

function stop(s::Pereira, to::TimerOutput, stats::AbstractStats, totalstats::AbstractStats)
    z_UB = stats.upperbound
    z_LB = stats.lowerbound
    K = stats.npaths
    σ = stats.σ_UB

    if totalstats.niterations > 0
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
