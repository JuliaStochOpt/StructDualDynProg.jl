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

function SOI.addcut!(sp::StochasticProgram, state, pool::SOI.AbstractSolutionPool, to::TimerOutput, ztol)
    if SOI.allfeasible(pool)
        gencut(SOI.get(sp, CutGenerator(), state), sp, state, pool, to::TimerOutput, ztol)
    else
        gencut(FeasibilityCutGenerator(), sp, state, pool, to::TimerOutput, ztol)
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
function SOI.addcut!(sp::StochasticProgram, state, cut::AveragedOptimalityCut)
    addcut(nodedata(sp, state).nlds.localOC, cut.a, cut.β, nodedata(sp, state).nlds)
end

function SOI.applycuts!(sp::StochasticProgram, state)
    applycut(FeasibilityCutGenerator(), sp, state)
    applycut(SOI.get(sp, CutGenerator(), state), sp, state)
end

function SOI.applycuts!(sp::StochasticProgram, tr::Transition, ::Type{<:FeasibilityCut})
    apply!(nodedata(sp, SOI.get(sp, SOI.Target(), tr)).fcuts)
end
function SOI.applycuts!(sp::StochasticProgram, tr::Transition, ::Type{<:MultiOptimalityCut})
    apply!(nodedata(sp, SOI.get(sp, SOI.Target(), tr)).ocuts)
end
function SOI.applycuts!(sp::StochasticProgram, state, ::Type{<:AveragedOptimalityCut})
    apply!(nodedata(sp, state).nlds.localOC)
end
