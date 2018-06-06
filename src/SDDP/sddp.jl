export SDDP

function solvejobs!(sp, jobsd, stats, stopatinf)
    infeasibility_detected = false
    for (node, jobs) in jobsd
        for job in jobs
            if !stopatinf || job.parent.pool.children_feasible
                solvejob!(sp, job, node, stats)
                if !job.parent.pool.children_feasible
                    infeasibility_detected = true
                end
            end
        end
    end
    infeasibility_detected
end

function gencuts(pathsd, sp, stats, ztol)
    for (parent, paths) in pathsd
        for path in paths
            if path.pool.children_feasible
                SOI.gencut(SOI.get(sp, SOI.CutGenerator(), parent), sp, parent, path.pool, stats, ztol)
            else
                SOI.gencut(SOI.FeasibilityCutGenerator(), sp, parent, path.pool, stats, ztol)
            end
        end
    end
end

function applycuts(pathsd, sp)
    for (node, _) in pathsd
        SOI.applycut(SOI.FeasibilityCutGenerator(), sp, node)
        SOI.applycut(SOI.get(sp, SOI.CutGenerator(), node), sp, node)
    end
end

"""
    struct Algorithm{ST<:AbstractPathSampler}
        K::Int
        verbose::Int
        pathsampler::ST
        ztol::Float64
        stopatinf::Bool
        mergepaths::Bool
        forwardcuts::Bool
        backwardcuts::Bool
    end

SDDP algorithm exploring `K` paths per iteration stages.
The paths will be selected according to `pathsampler` and equivalent paths might be merged if their difference is smaller than `ztol` and `mergepaths` is true.
The parameter `ztol` is also used to check whether a new cut is useful.
When a scenario is infeasible and `stopatinf` is true then no other scenario with the same ancestor is run. Note that since the order in which the different scenarios is run is not deterministic, this might introduce nondeterminism even if the sampling is deterministic.
By default, the cuts are added backward. However, if `forwardcuts` is set to `true` and `backwardcuts` is set to `false` the cuts are added forward.
"""
struct Algorithm{ST<:AbstractPathSampler} <: SOI.AbstractAlgorithm
    K::Int
    pathsampler::ST
    ztol::Float64
    stopatinf::Bool
    mergepaths::Bool
    forwardcuts::Bool
    backwardcuts::Bool
end
function Algorithm(; K::Int=25, pathsampler::AbstractPathSampler=ProbaPathSampler(true), ztol=1e-6, stopatinf=false, mergepaths=true, forwardcuts=false, backwardcuts=true)
    Algorithm(K, pathsampler, ztol, stopatinf, mergepaths, forwardcuts, backwardcuts)
end

function SOI.forwardpass!(sp::SOI.AbstractStochasticProgram, algo::Algorithm, verbose)
    stats = SOI.SDDPStats()

    master = SOI.get(sp, SOI.MasterState())
    NodeT = typeof(master)
    stats.solvertime += SOI.@_time mastersol = SOI.get(sp, SOI.Solution(), master)
    stats.nsolved += 1
    stats.niterations += 1
    infeasibility_detected = SOI.getstatus(mastersol) == :Infeasible

    num_stages = SOI.get(sp, SOI.NumberOfStages())
    PathT = SDDPPath{SOI.get(sp, SOI.TransitionType()), typeof(mastersol)}
    TT = Tuple{NodeT, Vector{PathT}}
    pathsd = Vector{Vector{TT}}(num_stages)
    if infeasibility_detected
        pathsd[1] = TT[]
    else
        pathsd[1] = TT[(master, [SDDPPath{SOI.get(sp, SOI.TransitionType())}(mastersol, [SOI.getstateobjectivevalue(mastersol)], [1.], [algo.K])])]
    end
    endedpaths = PathT[]

    for t in 2:num_stages
        verbose >= 3 && println("Stage $t/$num_stages")

        # Merge paths
        if algo.mergepaths
            pathsd[t-1] = mergesamepaths(pathsd[t-1], stats, algo.ztol)
        end

        # Make jobs
        jobsd = childjobs(sp, pathsd[t-1], algo.pathsampler, t, num_stages, endedpaths)

        # Solve Jobs (parallelism possible here)
        infeasibility_detected |= solvejobs!(sp, jobsd, stats, algo.stopatinf)

        if algo.forwardcuts
            gencuts(pathsd[t-1], sp, stats, algo.ztol)
            # The cut is added after so that they are added by group and duplicate detection in CutPruners works better
            applycuts(pathsd[t-1], sp)
        end

        # Jobs -> Paths
        pathsd[t] = jobstopaths(jobsd, sp)
    end

    if infeasibility_detected
        z_UB = Inf # FIXME assumes minimization
        σ = 0
    else
        for (_, paths) in pathsd[end]
            append!(endedpaths, paths)
        end
        z_UB, σ = meanstdpaths(endedpaths, algo.K)
    end

    # update stats
    stats.upperbound = z_UB
    stats.σ_UB = σ
    stats.npaths = algo.K
    stats.lowerbound = SOI.getobjectivevalue(mastersol)

    pathsd, mastersol, stats
end

function SOI.backwardpass!(sp::SOI.AbstractStochasticProgram, algo::Algorithm, pathsd, verbose)
    stats = SOI.SDDPStats()
    if algo.backwardcuts
        # We could ask a node if it has new fcuts/ocuts before resolving, that way the last stage won't be a particular case anymore
        num_stages = SOI.get(sp, SOI.NumberOfStages())
        for t in num_stages:-1:2
            # The last stages do not need to be resolved since it does not have any new cut
            if t != num_stages
                # Make jobs
                endedpaths = SDDPPath{SOI.get(sp, SOI.TransitionType())}[]
                jobsd = childjobs(sp, pathsd[t-1], algo.pathsampler, t, num_stages, endedpaths) # FIXME shouldn't need pathsampler here
                # Solve Jobs (parallelism possible here)
                solvejobs!(sp, jobsd, stats, algo.stopatinf)
            end

            gencuts(pathsd[t-1], sp, stats, algo.ztol)
            # The cut is added after so that they are added by group and duplicate detection in CutPruners works better
            applycuts(pathsd[t-1], sp)
        end
    end
    stats
end
