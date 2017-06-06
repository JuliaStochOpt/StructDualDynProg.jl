export model2lattice, getSDDPNode, SDDPclear

type SDDPModelData
    nodes::Vector{Nullable{SDDPNode}}
end

function getSDDPNode(m::Model, t, num_stages, solver, parent, pruningalgo::AbstractCutPruningAlgo, cutmode::Symbol, detectlb::Bool=true, newcut::Symbol=:InvalidateSolver)
    if !(:SDDP in keys(m.ext))
        nodes = Vector{Nullable{SDDPNode}}(num_stages)
        fill!(nodes, nothing)
        m.ext[:SDDP] = SDDPModelData(nodes)
    end
    nodes = m.ext[:SDDP].nodes
    if isnull(nodes[t])
        # The last argument contains the categories (e.g. :Cont, :Int, :Bool, ...) but it is currently unused
        c, T, W, h, C, K, _ = StructJuMP.conicconstraintdata(m)
        newnode = SDDPNode(NLDS{Float64}(W,h,T,K,C,c,solver,pruningalgo, newcut), parent)
        nodes[t] = newnode
        struc = getStructure(m)
        if t < num_stages
            num_scen = length(struc.children)
            children = Vector{SDDPNode{Float64}}(num_scen)
            probability = Vector{Float64}(num_scen)
            for (i, id) in enumerate(keys(struc.children))
                children[i] = getSDDPNode(struc.children[id], t+1, num_stages, solver, newnode, pruningalgo, cutmode, detectlb, newcut)
                probability[i] = struc.probability[id]
            end
            setchildren!(newnode, children, probability, cutmode)
            if detectlb
                setθlb!(newnode, map(getobjlb, children))
            end
        end
    end
    get(nodes[t])
end

"""
$(SIGNATURES)

Transforms a [StructJuMP](https://github.com/StructJuMP/StructJuMP.jl) model into a lattice that can be used by the SDDP algorithm.
The master problem is assumed to have model `m` and the scenarios are considered up to `num_stages` stages.
The `pruningalgo` is as defined in [CutPruners](https://github.com/JuliaPolyhedra/CutPruners.jl).
If `cutmode` is `:MultiCut`, one variable `θ_i` is created for each scenario. Otherwise, if `cutmode` is `:AveragedCut`, only one variable `θ` is created and it represents the expected value of the objective value of the scenarios. If `cutmode` is `:NoOptimalityCut` then no `θ` is created, only use this option if the objective of all models is zero except fo the master model.
"""
function model2lattice(m::Model, num_stages, solver, pruningalgo::AbstractCutPruningAlgo, cutmode::Symbol=:MultiCut, detectlb::Bool=true, newcut::Symbol=:InvalidateSolver)
    root = getSDDPNode(m, 1, num_stages, solver, nothing, pruningalgo, cutmode, detectlb, newcut)
	GraphSDDPTree(root)
end

function SDDPclear(m::Model)
    if :SDDP in keys(m.ext)
        pop!(m.ext, :SDDP)
        for (id, child) in getStructure(m).children
            SDDPclear(child)
        end
    end
end
