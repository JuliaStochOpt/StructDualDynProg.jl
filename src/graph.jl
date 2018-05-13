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

abstract type AbstractTransition end

# Mutable for setprobability!
mutable struct Transition{S} <: AbstractTransition
    source::Int
    target::Int
    σ::Int # FIXME NLDS is the only one who needs to map out_transitions to 1:n when it uses AveragedCut, it should have an internal dictionary
    proba::S
    childT::Nullable{AbstractMatrix{S}}
    function Transition(source::Int, target::Int, σ::Int, proba::S, childT) where S
        new{S}(source, target, σ, proba, childT)
    end
end
source(sp::AbstractStochasticProgram, tr::Transition) = tr.source
target(sp::AbstractStochasticProgram, tr::Transition) = tr.target
probability(sp::AbstractStochasticProgram, tr::Transition) = tr.proba

"""
    StochasticProgram{S, TT}

StochasticProgram of coefficient type `S` and transition type `TT`.
"""
mutable struct StochasticProgram{S, TT} <: AbstractStochasticProgram
    out_transitions::Vector{Vector{TT}} # out_transitions[i] : outgoing transitions for node i
    data::Vector{NodeData{S}}           # data[i] : data for node i
    function StochasticProgram{S, TT}() where {S, TT}
        new{S, TT}(Vector{TT}[], NodeData{S}[])
    end
end
StochasticProgram{S}() where S = StochasticProgram{S, Transition{S}}()

nodedata(sp::StochasticProgram, node::Int) = sp.data[node]

getobjectivebound(sp::StochasticProgram, node) = getobjectivebound(nodedata(sp, node).nlds)
setθbound!(sp::StochasticProgram, node, tr, θlb) = setθbound!(nodedata(sp, node).nlds, edgeid(sp, tr), θlb)
statedim(sp::StochasticProgram, node) = nodedata(sp, node).nlds.nx

# LightGraphs interface
out_transitions(sp::StochasticProgram, node::Int) = sp.out_transitions[node]
transitiontype(::StochasticProgram{S, TT}) where {S, TT} = TT
# May be different from the number of out-neighbors if there are multiple transitions with the same target
LightGraphs.outdegree(sp::StochasticProgram, node::Int) = length(out_transitions(sp, node))

getmaster(sp::StochasticProgram) = 1

# If the graph is not a tree, this will loop if I don't use a num_stages limit
function numberofpaths(sp::StochasticProgram, node, len)
    @assert len >= 0
    if iszero(len) || isleaf(sp, node)
        1
    else
        npath = nodedata(sp, node).npath
        if !(len in keys(npath))
            npath[len] = sum(map(tr -> numberofpaths(sp, target(sp, tr), len-1), out_transitions(sp, node)))
        end
        npath[len]
    end
end

cutgenerator(sp::StochasticProgram, node) = nodedata(sp, node).nlds.cutgen
function setcutgenerator!(sp::StochasticProgram, node, cutgen::AbstractOptimalityCutGenerator)
    nodedata(sp, node).nlds.cutgen = cutgen
end

function add_scenario_state!(sp::StochasticProgram{S}, data::NodeData) where S
    push!(sp.out_transitions, Transition{S}[])
    push!(sp.data, data)
    @assert length(sp.out_transitions) == length(sp.data)
    length(sp.data)
end

function add_scenario_transition!(sp::StochasticProgram, parent, child, proba, childT=nothing)
    tr = Transition(parent, child, outdegree(sp, parent)+1, proba, childT)
    push!(sp.out_transitions[parent], tr)
    data = nodedata(sp, parent)
    empty!(data.npath)
    childdata = nodedata(sp, child)
    add_scenario_transition!(data.nlds, childdata.fcuts, childdata.ocuts, proba, childT)
    @assert length(data.nlds.childFC) == length(data.nlds.proba) == outdegree(sp, parent)
    tr
end

function setprobability!(sp::StochasticProgram, tr, proba)
    tr.proba = proba
    data = nodedata(sp, source(sp, tr))
    setprobability!(data.nlds, edgeid(sp, tr), proba)
end

edgeid(sp::StochasticProgram, tr::Transition) = tr.σ

function solve!(sp::StochasticProgram, node)
    getsolution(nodedata(sp, node).nlds)
end

function setchildx!(sp::StochasticProgram, tr, sol::Solution)
    data = nodedata(sp, source(sp, tr))
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
    setparentx(nodedata(sp, target(sp, tr)).nlds, x, xuray, sol.objvalxuray)
end

function getθvalue(sp::StochasticProgram, tr::AbstractTransition, sol::Solution)
    @assert length(sol.θ) == outdegree(sp, source(sp, tr))
    getθvalue(sol, edgeid(sp, tr))
end

function getθvalue(sp::StochasticProgram, node, sol::Solution)
    @assert length(sol.θ) == 1
    getθvalue(sol, 1)
end

function add_feasibility_cut!(sp::StochasticProgram, node, coef, rhs, author)
    # coef is a ray
    # so alpha * coef is also valid for any alpha >= 0.
    # Hence coef might have very large coefficients and alter
    # the numerial accuracy of the master's solver.
    # We scale it to avoid this issue
    scaling = max(abs(rhs), maximum(abs, coef))
    addcut(nodedata(sp, node).fcuts, coef/scaling, sign(rhs), nodedata(sp, author).nlds)
end
function add_optimality_cut!(sp::StochasticProgram, node, coef, rhs, author)
    addcut(nodedata(sp, node).nlds.localOC, coef, rhs, nodedata(sp, author).nlds)
end
function add_optimality_cut_for_parent!(sp::StochasticProgram, node, coef, rhs, author)
    addcut(nodedata(sp, node).ocuts, coef, rhs, nodedata(sp, author).nlds)
end

function apply_feasibility_cuts!(sp::StochasticProgram, node)
    apply!(nodedata(sp, node).fcuts)
end
function apply_optimality_cuts!(sp::StochasticProgram, node)
    apply!(nodedata(sp, node).nlds.localOC)
end
function apply_optimality_cuts_for_parent!(sp::StochasticProgram, node)
    apply!(nodedata(sp, node).ocuts)
end
