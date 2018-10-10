module WaitAndSee

using Compat, Compat.SparseArrays

using StructDualDynProg
using StochOptInterface
const SOI = StochOptInterface

using MathProgBase

struct Algorithm{S<:MathProgBase.AbstractMathProgSolver} <: SOI.AbstractAlgorithm
    solver::S
    K::Int
end
Algorithm(solver::MathProgBase.AbstractMathProgSolver, K=25) = Algorithm(solver, K)

mutable struct WaitAndSeePath{NodeT}
    node::NodeT
    nlds::Vector{StructDualDynProg.StructProg.NLDS}
    z::Float64
    proba::Float64
    K::Int
end

function meanstdpaths(paths::Vector{WaitAndSeePath}, totalK)
    z = Float64[x.z for x in paths]
    proba = Float64[x.proba for x in paths]
    npaths = Int[x.K for x in paths]
    return StructDualDynProg.SDDP.meanstdpaths(z, proba, npaths, totalK)
end

function SOI.optimize!(sp::SOI.AbstractStochasticProgram, algo::Algorithm,
                       stopcrit::SOI.AbstractStoppingCriterion=SOI.IterLimit(0), # Not used
                       verbose::Int=0)
    master = SOI.get(sp, SOI.MasterNode())
    num_stages = SOI.get(sp, SOI.NumberOfStages())
    paths = WaitAndSeePath[WaitAndSeePath(master, StructDualDynProg.StructProg.NLDS[StructDualDynProg.StructProg.nodedata(sp, master).nlds], .0, 1., algo.K)]
    for t in 2:num_stages
        newpaths = WaitAndSeePath[]
        for path in paths
            if isempty(SOI.get(sp, SOI.OutTransitions(), path.node))
                push!(newpaths, path)
            else
                npaths = StructDualDynProg.SDDP.samplepaths(StructDualDynProg.SDDP.ProbaPathSampler(true), sp, path.node, path.K, t, num_stages)
                childs = algo.K == -1 ? (1:length(SOI.get(sp, SOI.OutTransitions(), path.node))) : findall(npaths .> 0)
                for (i, tr) in enumerate(SOI.get(sp, SOI.OutTransitions(), path.node))
                    if algo.K == -1 || npaths[i] > 0
                        push!(newpaths, WaitAndSeePath(SOI.get(sp, SOI.Target(), tr), [path.nlds; StructDualDynProg.StructProg.nodedata(sp, SOI.get(sp, SOI.Target(), tr)).nlds], path.z, path.proba * SOI.get(sp, SOI.Probability(), tr), npaths[i]))
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
        bs = Vector{Vector}(undef, length(path.nlds))
        Ks = []
        C = []
        c = similar(StructDualDynProg.StructProg.nodedata(sp, master).nlds.c, nvars[end])
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
            append!(C, [(cone, offsetnvars .+ idx) for (cone, idx) in nlds.C])
        end
        model = MathProgBase.LinearQuadraticModel(algo.solver)
        StructDualDynProg.StructProg._load!(model, c, A, bs, Ks, C)
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
    meanstdpaths(newpaths, algo.K)
end

end
