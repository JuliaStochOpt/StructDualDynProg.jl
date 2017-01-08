import Base.|, Base.&
export stop, AbstractStoppingCriterion, OrStoppingCriterion, AndStoppingCriterion, IterLimit, Pereira, CutLimit

abstract AbstractStoppingCriterion

"""
`stop(stopcrit, iter, nfcuts, nocuts, K, z_LB, z_UB, σ)`

Returns whether the SDDP algorithm should stop.
If `iter` is 0, no iteration has already been done, otherwise, the `iter`th iteration has just finished.
This iteration used `K` paths and generated `nfcuts` (resp. `nocuts`) new feasibility (resp. optimality) cuts.
The lower bound is now `z_LB` and the upper bound has mean `z_UB` and variance `σ`.
"""
function stop(s::AbstractStoppingCriterion, iter, nfcuts, nocuts, K, z_LB, z_UB, σ)
    error("`stop' function not defined for $(typeof(s))")
end

"""
$(TYPEDEF)

Stops if `lhs` *or* `rhs` want to stop.
"""
type OrStoppingCriterion <: AbstractStoppingCriterion
    lhs::AbstractStoppingCriterion
    rhs::AbstractStoppingCriterion
end

function stop(s::OrStoppingCriterion, iter, nfcuts, nocuts, K, z_LB, z_UB, σ)
    stop(s.lhs, iter, nfcuts, nocuts, K, z_LB, z_UB, σ) || stop(s.rhs, iter, nfcuts, nocuts, K, z_LB, z_UB, σ)
end

function (|)(lhs::AbstractStoppingCriterion, rhs::AbstractStoppingCriterion)
    OrStoppingCriterion(lhs, rhs)
end

"""
$(TYPEDEF)

Stops if `lhs` *and* `rhs` want to stop.
"""
type AndStoppingCriterion <: AbstractStoppingCriterion
    lhs::AbstractStoppingCriterion
    rhs::AbstractStoppingCriterion
end

function stop(s::AndStoppingCriterion, iter, nfcuts, nocuts, K, z_LB, z_UB, σ)
    stop(s.lhs, iter, nfcuts, nocuts, K, z_LB, z_UB, σ) && stop(s.rhs, iter, nfcuts, nocuts, K, z_LB, z_UB, σ)
end

function (&)(lhs::AbstractStoppingCriterion, rhs::AbstractStoppingCriterion)
    AndStoppingCriterion(lhs, rhs)
end

"""
$(TYPEDEF)

Stops if `iter` ≧ `limit`.
"""
type IterLimit <: AbstractStoppingCriterion
    limit::Int
end

function stop(s::IterLimit, iter, nfcuts, nocuts, K, z_LB, z_UB, σ)
    iter >= s.limit
end

"""
$(TYPEDEF)

Stops if there was less than or equal to `limit` cuts added in the iteration.
For instance, `CutLimit(0)` stops when there are no cuts added.
"""
type CutLimit <: AbstractStoppingCriterion
    limit::Int
end

function stop(s::CutLimit, iter, nfcuts, nocuts, K, z_LB, z_UB, σ)
    iter > 0 && nfcuts + nocuts <= s.limit
end

"""
$(TYPEDEF)

Stops if `z_UB - α * σ/√K < z_LB < z_UB - α * σ/√K` and `σ / √K > β * max(1, |z_LB|))`
"""
type Pereira <: AbstractStoppingCriterion
    α::Float64
    β::Float64

    Pereira(α=2.0, β=0.05) = new(α, β)
end

function stop(s::Pereira, iter, nfcuts, nocuts, K, z_LB, z_UB, σ)
    if iter > 0
        @assert K >= 0
        σ1 = σ / √K
        σ2 = s.α * σ1
        z_UB - σ2 <= z_LB <= z_UB + σ2 && σ1 < s.β * max(1, abs(z_LB))
    else
        false
    end
end
