export AbstractSDDPTree
@compat abstract type AbstractSDDPTree{S} end

type GraphSDDPTree{S} <: AbstractSDDPTree{S}
	root::SDDPNode{S}
end

getmaster(g::GraphSDDPTree) = g.root, g.root

# Get children scenarios
haschildren(g::GraphSDDPTree, node::SDDPNode) = !isempty(node.children)
nchildren(g::GraphSDDPTree, node::SDDPNode) = length(node.children)
children(g::GraphSDDPTree, node::SDDPNode) = node.children
getchild(g::GraphSDDPTree, node::SDDPNode, i) = node.children[i]
# Get proba of children scenario
getproba(g::GraphSDDPTree, node::SDDPNode, i) = node.proba[i]
getprobas(g::GraphSDDPTree, node::SDDPNode) = node.proba

# Get number of paths
numberofpaths(g::GraphSDDPTree, node::SDDPNode, t, num_stages) = numberofpaths(node, t, num_stages)

cutmode(g::GraphSDDPTree, node::SDDPNode) = node.nlds.cutmode
