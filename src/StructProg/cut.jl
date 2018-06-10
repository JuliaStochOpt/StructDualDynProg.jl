abstract type AbstractCut end

"""
    addcut!(sp::AbstractStochasticProgram, state_or_tr, cut::AbstractCut)

Add cut `cut` to the state or transition `state_or_tr`.
"""
function addcut! end

"""
    applycuts!(sp::AbstractStochasticProgram, state, T::Type{<:AbstractCut})

Apply cuts of type `T` to the state or transition `state_or_tr`.
"""
function applycuts! end

#"""
#    iscut(sol::AbstractSolution, cut::AbstractCut)
#
#Return a `Bool` indicating whether the solution `sol` is cut by the cut `cut`.
#That is, whether `sol` would still be a feasible solution with the same objective value if `cut` was added.
#"""
#function iscut end

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
