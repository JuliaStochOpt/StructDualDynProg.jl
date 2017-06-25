export AbstractSDDPGraph, haschildren, nchidren, children, getchild, getproba, getprobas, cutgen, numberofpaths
abstract type AbstractSDDPGraph{S} end

mutable struct SDDPGraph{S} <: AbstractSDDPGraph{S}
	root::SDDPNode{S}
end

getmaster(g::SDDPGraph) = g.root, g.root

# Get children scenarios
haschildren(g::SDDPGraph, node::SDDPNode) = !isempty(node.children)
nchildren(g::SDDPGraph, node::SDDPNode) = length(node.children)
children(g::SDDPGraph, node::SDDPNode) = node.children
getchild(g::SDDPGraph, node::SDDPNode, i) = node.children[i]
# Get proba of children scenario
getproba(g::SDDPGraph, node::SDDPNode, i) = node.proba[i]
getprobas(g::SDDPGraph, node::SDDPNode) = node.proba

# Get number of paths
numberofpaths(g::SDDPGraph, num_stages) = numberofpaths(g.root, 1, num_stages)
numberofpaths(g::SDDPGraph, node::SDDPNode, t, num_stages) = numberofpaths(node, t, num_stages)

cutgen(g::SDDPGraph, node::SDDPNode) = node.nlds.cutgen
