export AbstractSDDPTree
abstract AbstractSDDPTree{S}

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

cutmode(g::GraphSDDPTree, node::SDDPNode) = node.nlds.cutmode

function choosepaths(g::GraphSDDPTree, node::SDDPNode, K::Int, pathsel, t, num_stages)
    if K == -1
        map(child->-1, node.children)
    else
        if pathsel == :nPaths
            den = numberofpaths(node, t-1, num_stages)
            pmf = map(child->numberofpaths(child, t, num_stages) / den, node.children)
        else
            pmf = node.proba
        end
        cmf = cumsum(pmf)
        @assert abs(cmf[end] - 1) < 1e-6
        cmf[end] = 1
        samples = rand(Float64, K)
        npaths = zeros(Int, length(node.children))
        sort!(samples)
        i = 1
        for j in samples
            while j >= cmf[i]
                i += 1
            end
            npaths[i] += 1
        end
        npaths
    end
end

function choosepaths(g::GraphSDDPTree, node::SDDPNode, K::Vector{Int}, pathsel, t, num_stages)
    npathss = Vector{Int}[similar(K) for i in 1:length(node.children)]
    for i in 1:length(K)
        npaths = choosepaths(g, node, K[i], pathsel, t, num_stages)
        for c in 1:length(node.children)
            npathss[c][i] = npaths[c]
        end
    end
    npathss
end
