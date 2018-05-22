export waitandsee

mutable struct WaitAndSeePath{NodeT}
    node::NodeT
    nlds::Vector{NLDS}
    z::Float64
    proba::Float64
    K::Int
end

function meanstdpaths(paths::Vector{WaitAndSeePath}, totalK)
    z = Float64[x.z for x in paths]
    proba = Float64[x.proba for x in paths]
    npaths = Int[x.K for x in paths]
    meanstdpaths(z, proba, npaths, totalK)
end

function waitandsee(sp::SOI.AbstractStochasticProgram, num_stages, solver, totalK=25, verbose=0)
    master = SOI.get(sp, SOI.MasterState())
    paths = WaitAndSeePath[WaitAndSeePath(master, NLDS[nodedata(sp, master).nlds], .0, 1., totalK)]
    for t in 2:num_stages
        newpaths = WaitAndSeePath[]
        for path in paths
            if iszero(outdegree(sp, path.node))
                push!(newpaths, path)
            else
                npaths = samplepaths(ProbaPathSampler(true), sp, path.node, path.K, t, num_stages)
                childs = totalK == -1 ? (1:outdegree(sp, path.node)) : find(npaths .> 0)
                for (i, tr) in enumerate(SOI.get(sp, SOI.OutTransitions(), path.node))
                    if totalK == -1 || npaths[i] > 0
                        push!(newpaths, WaitAndSeePath(SOI.target(sp, tr), [path.nlds; nodedata(sp, SOI.target(sp, tr)).nlds], path.z, path.proba * SOI.get(sp, SOI.Probability(), tr), npaths[i]))
                    end
                end
            end
        end
        paths = newpaths
    end

    sumz = 0
    newpaths = WaitAndSeePath[]
    for path in paths
        nvars = cumsum(Int[nlds.nx for nlds in path.nlds])
        ncons = cumsum(Int[nlds.nπ for nlds in path.nlds])
        # b - Ax in K_1
        # x \in K_2
        A = spzeros(ncons[end], nvars[end])
        bs = Vector{Vector}(length(path.nlds))
        Ks = []
        C = []
        c = similar(nodedata(sp, master).nlds.c, nvars[end])
        for i in 1:length(path.nlds)
            offsetnvars = i == 1 ? 0 : nvars[i-1]
            offsetncons = i == 1 ? 0 : ncons[i-1]
            nlds = path.nlds[i]
            c[offsetnvars+1:nvars[i]] = nlds.c
            bs[i] = nlds.h
            A[offsetncons+1:ncons[i], offsetnvars+1:nvars[i]] = nlds.W
            if i > 1
                offset = i == 2 ? 0 : nvars[i-2]
                A[offsetncons+1:ncons[i], offset+1:nvars[i-1]] = nlds.T
            end
            push!(Ks, nlds.K)
            append!(C, [(cone, offsetnvars+idx) for (cone, idx) in nlds.C])
        end
        model = MathProgBase.LinearQuadraticModel(solver)
        _load!(model, c, A, bs, Ks, C)
        MathProgBase.optimize!(model)
        status = MathProgBase.status(model)
        if status == :Error
            error("The solver reported an error")
        elseif status == :UserLimit
            error("The solver reached iteration limit or timed out")
        elseif status == :Infeasible
            error("The problem is infeasible for some realization of the uncertainty")
        elseif status == :Unbounded
            error("The problem is unbounded for some realization of the uncertainty")
        else
            @assert status == :Optimal
            objval = MathProgBase.getobjval(model)
            push!(newpaths, WaitAndSeePath(path.node, path.nlds, objval, path.proba, path.K))
        end
    end
    meanstdpaths(newpaths, totalK)
end
