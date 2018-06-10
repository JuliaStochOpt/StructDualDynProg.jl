using LightGraphs

import MathProgBase

mutable struct NodeData{S}
    nlds::NLDS{S}
    npath::Dict{Int, Int}

    # Feasibility cuts
    fcuts::CutStore{S}
    # Optimality cuts
    ocuts::CutStore{S}

    function NodeData{S}(nlds::NLDS{S}, nvars_a) where S
        new{S}(nlds, Dict{Int, Int}(), CutStore{S}(nvars_a), CutStore{S}(nvars_a))
    end
end

NodeData(nlds::NLDS{S}, parent) where {S} = NodeData{S}(nlds, parent)

function Base.show(io::IO, data::NodeData)
    println(io, "Node of $(data.nlds.nx) variables")
end

# Mutable for setprobability!
mutable struct Transition{S} <: SOI.AbstractTransition
    source::Int
    target::Int
    σ::Int # FIXME NLDS is the only one who needs to map out_transitions to 1:n when it uses AveragedCut, it should have an internal dictionary
    proba::S
    childT::Nullable{AbstractMatrix{S}}
    function Transition(source::Int, target::Int, σ::Int, proba::S, childT) where S
        new{S}(source, target, σ, proba, childT)
    end
end
SOI.get(::SOI.AbstractStochasticProgram, ::SOI.Source, tr::Transition) = tr.source
SOI.get(::SOI.AbstractStochasticProgram, ::SOI.Target, tr::Transition) = tr.target
SOI.get(::SOI.AbstractStochasticProgram, ::SOI.Probability, tr::Transition) = tr.proba

Base.:(==)(t1::Transition, t2::Transition) = t1.source == t2.source && t1.target == t2.target && t1.σ == t2.σ
Base.hash(t::Transition, u::UInt) = hash(t.source, hash(t.target, hash(t.σ, u)))

"""
    StochasticProgram{S, TT}

StochasticProgram of coefficient type `S` and transition type `TT`.
"""
mutable struct StochasticProgram{S, TT} <: SOI.AbstractStochasticProgram
    out_transitions::Vector{Vector{TT}} # out_transitions[i] : outgoing transitions for node i
    data::Vector{NodeData{S}}           # data[i] : data for node i
    num_stages::Int
    function StochasticProgram{S, TT}(num_stages) where {S, TT}
        new{S, TT}(Vector{TT}[], NodeData{S}[], num_stages)
    end
end
StochasticProgram{S}(num_stages) where S = StochasticProgram{S, Transition{S}}(num_stages)

nodedata(sp::StochasticProgram, node::Int) = sp.data[node]

SOI.get(sp::StochasticProgram, ::SOI.NumberOfStages) = sp.num_stages

SOI.get(sp::StochasticProgram, ::SOI.StateObjectiveValueBound, state) = getobjectivebound(nodedata(sp, state).nlds)
SOI.set!(sp::StochasticProgram, ::SOI.TransitionObjectiveValueBound, tr::Transition, θlb) = setθbound!(nodedata(sp, SOI.get(sp, SOI.Source(), tr)).nlds, edgeid(sp, tr), θlb)
SOI.get(sp::StochasticProgram, ::SOI.Dimension, state) = nodedata(sp, state).nlds.nx

SOI.get(sp::StochasticProgram, ::SOI.OutTransitions, node::Int) = sp.out_transitions[node]
SOI.get(::StochasticProgram{S, TT}, ::SOI.TransitionType) where {S, TT} = TT

SOI.get(::StochasticProgram, ::SOI.MasterState) = 1

# If the graph is not a tree, this will loop if I don't use a num_stages limit
function SOI.get(sp::StochasticProgram, nop::SOI.NumberOfPathsFrom, node)
    len = nop.length
    @assert len >= 0
    if iszero(len) || iszero(outdegree(sp, node))
        1
    else
        npath = nodedata(sp, node).npath
        if !(len in keys(npath))
            npath[len] = sum(map(tr -> SOI.get(sp, SOI.NumberOfPathsFrom(len-1), SOI.get(sp, SOI.Target(), tr)), SOI.get(sp, SOI.OutTransitions(), node)))
        end
        npath[len]
    end
end

SOI.get(sp::StochasticProgram, ::SOI.CutGenerator, node) = nodedata(sp, node).nlds.cutgen
function SOI.set!(sp::StochasticProgram, ::SOI.CutGenerator, node, cutgen::SOI.AbstractOptimalityCutGenerator)
    nodedata(sp, node).nlds.cutgen = cutgen
end

function SOI.add_scenario_state!(sp::StochasticProgram{S}, data::NodeData) where S
    push!(sp.out_transitions, Transition{S}[])
    push!(sp.data, data)
    @assert length(sp.out_transitions) == length(sp.data)
    length(sp.data)
end

function SOI.add_scenario_transition!(sp::StochasticProgram, parent, child, proba, childT=nothing)
    tr = Transition(parent, child, outdegree(sp, parent)+1, proba, childT)
    push!(sp.out_transitions[parent], tr)
    data = nodedata(sp, parent)
    empty!(data.npath)
    childdata = nodedata(sp, child)
    add_scenario_transition!(data.nlds, childdata.fcuts, childdata.ocuts, proba, childT)
    @assert length(data.nlds.childFC) == length(data.nlds.proba) == outdegree(sp, parent)
    tr
end

function SOI.set!(sp::StochasticProgram, ::SOI.Probability, tr, proba)
    tr.proba = proba
    data = nodedata(sp, SOI.get(sp, SOI.Source(), tr))
    setprobability!(data.nlds, edgeid(sp, tr), proba)
end

edgeid(sp::StochasticProgram, tr::Transition) = tr.σ

function SOI.get(sp::StochasticProgram, ::SOI.Solution, node)
    getsolution(nodedata(sp, node).nlds)
end

function SOI.set!(sp::StochasticProgram, ::SOI.SourceSolution, tr, sol::Solution)
    data = nodedata(sp, SOI.get(sp, SOI.Source(), tr))
    if !isnull(tr.childT)
        T = get(tr.childT)
        x = T * sol.x
        xuray = sol.xuray
        if xuray !== nothing
            xuray = T * sol.xuray
        end
    else
        x = sol.x
        xuray = sol.xuray
    end
    setparentx(nodedata(sp, SOI.get(sp, SOI.Target(), tr)).nlds, x, xuray, sol.objvalxuray)
end

function SOI.getθvalue(sp::StochasticProgram, tr::SOI.AbstractTransition, sol::Solution)
    @assert length(sol.θ) == outdegree(sp, SOI.get(sp, SOI.Source(), tr))
    SOI.getθvalue(sol, edgeid(sp, tr))
end

function SOI.getθvalue(sp::StochasticProgram, node, sol::Solution)
    @assert length(sol.θ) == 1
    SOI.getθvalue(sol, 1)
end

function SOI.addcut!(sp::StochasticProgram, tr::Transition, cut::SOI.FeasibilityCut)
    # coef is a ray
    # so alpha * coef is also valid for any alpha >= 0.
    # Hence coef might have very large coefficients and alter
    # the numerial accuracy of the master's solver.
    # We scale it to avoid this issue
    scaling = max(abs(cut.β), maximum(abs, cut.a))
    addcut(nodedata(sp, SOI.get(sp, SOI.Target(), tr)).fcuts, cut.a/scaling, sign(cut.β), nodedata(sp, SOI.get(sp, SOI.Source(), tr)).nlds)
end
function SOI.addcut!(sp::StochasticProgram, tr::Transition, cut::SOI.MultiOptimalityCut)
    addcut(nodedata(sp, SOI.get(sp, SOI.Target(), tr)).ocuts, cut.a, cut.β, nodedata(sp, SOI.get(sp, SOI.Source(), tr)).nlds)
end
function SOI.addcut!(sp::StochasticProgram, state, cut::SOI.AveragedOptimalityCut)
    addcut(nodedata(sp, state).nlds.localOC, cut.a, cut.β, nodedata(sp, state).nlds)
end

function SOI.applycuts!(sp::StochasticProgram, tr::Transition, ::Type{<:SOI.FeasibilityCut})
    apply!(nodedata(sp, SOI.get(sp, SOI.Target(), tr)).fcuts)
end
function SOI.applycuts!(sp::StochasticProgram, tr::Transition, ::Type{<:SOI.MultiOptimalityCut})
    apply!(nodedata(sp, SOI.get(sp, SOI.Target(), tr)).ocuts)
end
function SOI.applycuts!(sp::StochasticProgram, state, ::Type{<:SOI.AveragedOptimalityCut})
    apply!(nodedata(sp, state).nlds.localOC)
end
