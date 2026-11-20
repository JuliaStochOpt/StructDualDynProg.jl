module WaitAndSee

using JuMP

using StructDualDynProg
using StochOptInterface
const SOI = StochOptInterface

struct Algorithm{S} <: SOI.AbstractAlgorithm
    solver::S
    K::Int
end
Algorithm(solver) = Algorithm(solver, 25)

mutable struct WaitAndSeePath{NodeT}
    node::NodeT
    transitions::Vector{StructDualDynProg.StructProg.Transition}
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
    paths = WaitAndSeePath[WaitAndSeePath(master, StructDualDynProg.StructProg.Transition[], .0, 1., algo.K)]
    for t in 2:num_stages
        newpaths = WaitAndSeePath[]
        for path in paths
            if isempty(SOI.get(sp, SOI.OutTransitions(), path.node))
                push!(newpaths, path)
            else
                npaths = StructDualDynProg.SDDP.samplepaths(StructDualDynProg.SDDP.ProbaPathSampler(true), sp, path.node, path.K, t, num_stages)
                for (i, tr) in enumerate(SOI.get(sp, SOI.OutTransitions(), path.node))
                    if algo.K == -1 || npaths[i] > 0
                        push!(newpaths, WaitAndSeePath(SOI.get(sp, SOI.Target(), tr), [path.transitions; tr], path.z, path.proba * SOI.get(sp, SOI.Probability(), tr), npaths[i]))
                    end
                end
            end
        end
        paths = newpaths
    end

    newpaths = WaitAndSeePath[]
    for path in paths
        model = JuMP.Model(algo.solver)
        root = StructDualDynProg.StructProg.nodedata(sp, master).nlds
        variables, objective = StructDualDynProg.StructProg.add_stage!(model, root.stage, JuMP.VariableRef[])
        for tr in path.transitions
            parent = isnothing(tr.childT) ? variables : tr.childT * variables
            child = StructDualDynProg.StructProg.nodedata(sp, tr.target).nlds
            variables, cost = StructDualDynProg.StructProg.add_stage!(model, child.stage, parent)
            JuMP.add_to_expression!(objective, cost)
        end
        @objective(model, Min, objective)
        JuMP.optimize!(model)
        status = JuMP.termination_status(model)
        if status == MOI.OPTIMAL
            objval = JuMP.objective_value(model)
            push!(newpaths, WaitAndSeePath(path.node, path.transitions, objval, path.proba, path.K))
        else
            error("The solver terminated with $status")
        end
    end
    meanstdpaths(newpaths, algo.K)
end

end
