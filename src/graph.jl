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

ET = LightGraphs.SimpleGraphs.SimpleEdge{Int64}

mutable struct StochasticProgram{S} <: AbstractStochasticProgram
    graph::LightGraphs.SimpleGraphs.SimpleDiGraph{Int}
    eid::Dict{ET, Int}
    proba::Dict{ET, S}
    childT::Dict{ET, AbstractMatrix{S}}
    data::Vector{NodeData{S}}
    function StochasticProgram{S}() where S
        new{S}(DiGraph(), Dict{ET, Int}(), Dict{ET, S}(), Dict{ET, AbstractMatrix{S}}(), NodeData{S}[])
    end
end
nodedata(sp::StochasticProgram, node::Int) = sp.data[node]

getobjlb(sp::StochasticProgram, node) = getobjlb(nodedata(sp, node).nlds)
setθlb!(sp::StochasticProgram, node, θlb) = setθlb!(nodedata(sp, node).nlds, θlb)
statedim(sp::StochasticProgram, node) = nodedata(sp, node).nlds.nx

# LightGraphs interface
LightGraphs.out_neighbors(sp::StochasticProgram, node::Int) = out_neighbors(sp.graph, node)

getmaster(sp::StochasticProgram) = 1

# If the graph is not a tree, this will loop if I don't use a num_stages limit
function numberofpaths(sp::StochasticProgram, node, len)
    @assert len >= 0
    if iszero(len) || isleaf(sp, node)
        1
    else
        npath = nodedata(sp, node).npath
        if !(len in keys(npath))
            npath[len] = sum(map(c -> numberofpaths(sp, c, len-1), out_neighbors(sp, node)))
        end
        npath[len]
    end
end

cutgenerator(sp::StochasticProgram, node) = nodedata(sp, node).nlds.cutgen
function setcutgenerator!(sp::StochasticProgram, node, cutgen::AbstractOptimalityCutGenerator)
    nodedata(sp, node).nlds.cutgen = cutgen
end

function add_scenario_state!(sp::StochasticProgram, data::NodeData)
    @assert add_vertex!(sp.graph)
    push!(sp.data, data)
    @assert nv(sp.graph) == length(sp.data)
    length(sp.data)
end

function add_scenario_transition!(sp::StochasticProgram, parent, child, proba, childT=nothing)
    edge = Edge(parent, child)
    if !add_edge!(sp.graph, edge)
        error("Edge already in the graph, multiple edges not supported yet")
    end
    @assert !haskey(sp.eid, edge)
    @assert !haskey(sp.proba, edge)
    @assert !haskey(sp.childT, edge)
    data = nodedata(sp, parent)
    sp.eid[edge] = outdegree(sp, parent)
    sp.proba[edge] = proba
    if childT !== nothing
        sp.childT[edge] = childT
    end
    empty!(data.npath)
    childdata = nodedata(sp, child)
    add_scenario_transition!(data.nlds, childdata.fcuts, childdata.ocuts, proba, childT)
    @assert length(data.nlds.childFC) == length(data.nlds.proba) == outdegree(sp, parent)
end

probability(sp::StochasticProgram, edge) = sp.proba[edge]

function setprobability!(sp::StochasticProgram, edge, proba)
    sp.proba[edge] = proba
    data = nodedata(sp, src(edge))
    setprobability!(data.nlds, edgeid(sp, edge), proba)
end

function edgeid(sp::StochasticProgram, edge)
    sp.eid[edge]
end

function solve!(sp::StochasticProgram, node)
    getsolution(nodedata(sp, node).nlds)
end

function setchildx!(sp::StochasticProgram, node, child, sol::Solution)
    data = nodedata(sp, node)
    edge = Edge(node, child)
    if haskey(sp.childT, edge)
        T = data.childT[edge]
        x = T * sol.x
        if sol.xuray !== nothing
            xuray = T * sol.xuray
        end
    else
        x = sol.x
        xuray = sol.xuray
    end
    setparentx(nodedata(sp, child).nlds, x, xuray, sol.objvalxuray)
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
