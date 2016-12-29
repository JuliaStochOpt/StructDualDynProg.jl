import MathProgBase

export SDDPNode, setchildren!, appendchildren!

type SDDPNode{S}
    nlds::NLDS{S}
    nvars::Int
    parent::Nullable{SDDPNode{S}}
    children::Vector{SDDPNode{S}}
    proba::Vector{S}
    childT::Nullable{Vector{AbstractMatrix{S}}}
    root::Bool
    leaf::Bool

    npath::Dict{Tuple{Int,Int},Int}

    # Feasibility cuts
    fcuts::CutStore{S}
    # Optimality cuts
    ocuts::CutStore{S}

    function SDDPNode(nlds::NLDS{S}, parent)
        nvars = size(nlds.W, 2)
        root = parent === nothing
        nvars_a = root ? 0 : parent.nvars
        new(nlds, nvars, parent, SDDPNode[], Float64[], nothing, root, true, Dict{Tuple{Int,Int},Int}(), CutStore{S}(nvars_a), CutStore{S}(nvars_a))
    end

end

SDDPNode{S}(nlds::NLDS{S}, parent) = SDDPNode{S}(nlds, parent)

function Base.show(io::IO, node::SDDPNode)
    if node.root
        print(io, "Root ")
        if node.leaf
            print(io, "and leaf ")
        end
        print(io, "n")
    else
        print(io, "N")
    end
    println(io, "ode of $(node.nvars) variables and outdegree of $(length(node.children)) with proba:")
    println(io, node.proba)
end

function setchildren!(node::SDDPNode, children, proba, cutmode, childT=nothing)
    @assert length(children) == length(proba)
    node.children = children
    node.proba = proba
    node.leaf = false
    childFC = map(child -> child.fcuts, children)
    childOC = map(child -> child.ocuts, children)
    node.childT = childT
    empty!(node.npath)
    setchildren!(node.nlds, childFC, childOC, proba, cutmode, childT)
end

function appendchildren!(node::SDDPNode, children, proba, childT=nothing)
    append!(node.children, children)
    if length(proba) == length(children)
        append!(node.proba, proba)
    else
        @assert length(proba) == length(node.children)
        node.proba = proba
    end
    node.leaf = false
    childFC = map(child -> child.fcuts, children)
    childOC = map(child -> child.ocuts, children)
    if childT === nothing
        @assert isnull(node.childT)
    else
        # If there isn't any child yet, node.childT is null
        if isnull(node.childT)
            node.childT = childT
        else
            append!(get(node.childT), childT)
        end
    end
    empty!(node.npath)
    appendchildren!(node.nlds, childFC, childOC, node.proba, childT)
    @assert length(node.nlds.childFC) == length(node.children)
end

function setchildx(node::SDDPNode, i::Int, sol::NLDSSolution)
    if !isnull(node.childT)
        x = get(node.childT)[i] * sol.x
        if sol.xuray !== nothing
            xuray = get(node.childT)[i] * sol.xuray
        end
    else
        x = sol.x
        xuray = sol.xuray
    end
    setparentx(node.children[i].nlds, x, xuray, sol.objvalxuray)
end

# If the graph is not a tree, this will loop if I don't use a num_stages limit
function numberofpaths(node::SDDPNode, t, num_stages)
    if t == num_stages || isempty(node.children)
        1
    else
        key = (t,num_stages)
        if !(key in keys(node.npath))
            node.npath[key] = sum(map(c -> numberofpaths(c, t+1, num_stages), node.children))
        end
        node.npath[key]
    end
end

function pushfeasibilitycut!(node, coef, rhs, author)
    # coef is a ray
    # so alpha * coef is also valid for any alpha >= 0.
    # Hence coef might have very large coefficients and alter
    # the numerial accuracy of the master's solver.
    # We scale it to avoid this issue
    scaling = max(abs(rhs), maximum(abs(coef)))
    addcut(node.fcuts, coef/scaling, sign(rhs), author.nlds)
end
function pushoptimalitycut!(node, coef, rhs, author)
    addcut(node.nlds.localOC, coef, rhs, author.nlds)
end
function pushoptimalitycutforparent!(node, coef, rhs, author)
    addcut(node.ocuts, coef, rhs, author.nlds)
end

function applyfeasibilitycut!(node)
    apply!(node.fcuts)
end
function applyoptimalitycut!(node)
    apply!(node.nlds.localOC)
end
function applyoptimalitycutforparent!(node)
    apply!(node.ocuts)
end

function loadAndSolve(node::SDDPNode)
    getsolution(node.nlds)
end
