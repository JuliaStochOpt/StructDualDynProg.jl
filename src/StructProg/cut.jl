abstract type AbstractCut end

struct FeasibilityCut{T, VT<:AbstractVector{T}} <: AbstractCut
    a::VT
    β::T
end

abstract type AbstractOptimalityCut <: AbstractCut end

struct MultiOptimalityCut{T, VT<:AbstractVector{T}} <: AbstractCut
    a::VT
    β::T
end

struct AveragedOptimalityCut{T, VT<:AbstractVector{T}} <: AbstractCut
    a::VT
    β::T
end

function SOI.addcut!(sp::StochasticProgram, node, pool::SOI.AbstractSolutionPool, to::TimerOutput, ztol)
    if SOI.allfeasible(pool)
        gencut(SOI.get(sp, CutGenerator(), node), sp, node, pool, to::TimerOutput, ztol)
    else
        gencut(FeasibilityCutGenerator(), sp, node, pool, to::TimerOutput, ztol)
    end
end

function SOI.addcut!(sp::StochasticProgram, tr::Transition, cut::FeasibilityCut)
    # coef is a ray
    # so alpha * coef is also valid for any alpha >= 0.
    # Hence coef might have very large coefficients and alter
    # the numerial accuracy of the master's solver.
    # We scale it to avoid this issue
    scaling = max(abs(cut.β), maximum(abs, cut.a))
    addcut(nodedata(sp, SOI.get(sp, SOI.Target(), tr)).fcuts, cut.a/scaling, sign(cut.β), nodedata(sp, SOI.get(sp, SOI.Source(), tr)).nlds)
end
function SOI.addcut!(sp::StochasticProgram, tr::Transition, cut::MultiOptimalityCut)
    addcut(nodedata(sp, SOI.get(sp, SOI.Target(), tr)).ocuts, cut.a, cut.β, nodedata(sp, SOI.get(sp, SOI.Source(), tr)).nlds)
end
function SOI.addcut!(sp::StochasticProgram, node, cut::AveragedOptimalityCut)
    addcut(nodedata(sp, node).nlds.localOC, cut.a, cut.β, nodedata(sp, node).nlds)
end

function SOI.applycuts!(sp::StochasticProgram, node)
    applycut(FeasibilityCutGenerator(), sp, node)
    applycut(SOI.get(sp, CutGenerator(), node), sp, node)
end

function SOI.applycuts!(sp::StochasticProgram, tr::Transition, ::Type{<:FeasibilityCut})
    apply!(nodedata(sp, SOI.get(sp, SOI.Target(), tr)).fcuts)
end
function SOI.applycuts!(sp::StochasticProgram, tr::Transition, ::Type{<:MultiOptimalityCut})
    apply!(nodedata(sp, SOI.get(sp, SOI.Target(), tr)).ocuts)
end
function SOI.applycuts!(sp::StochasticProgram, node, ::Type{<:AveragedOptimalityCut})
    apply!(nodedata(sp, node).nlds.localOC)
end
