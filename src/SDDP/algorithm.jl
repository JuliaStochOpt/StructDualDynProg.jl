function solvejobs!(sp, jobsd, to::TimerOutput, stopatinf::Bool)
    infeasibility_detected = false
    for (node, jobs) in jobsd
        for job in jobs
            if !stopatinf || SOI.allfeasible(job.parent.pool)
                solvejob!(sp, job, node, to)
                if !SOI.allfeasible(job.parent.pool)
                    infeasibility_detected = true
                end
            end
        end
    end
    infeasibility_detected
end

function gencuts(pathsd, sp, to::TimerOutput, ztol)
    for (parent, paths) in pathsd
        for path in paths
            SOI.addcut!(sp, parent, path.pool, to, ztol)
        end
    end
end

function applycuts(pathsd, sp)
    for (node, _) in pathsd
        SOI.applycuts!(sp, node)
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

struct Paths{NodeT, TT, SolT} <: SOI.AbstractPaths
    paths::Vector{Vector{Tuple{NodeT, Vector{SDDPPath{TT, SolT}}}}}
end
function SOI.npaths(paths::Paths)
    # At the master node, all paths are still concentrated in the same one
    node_paths = paths.paths[1]
    @assert length(node_paths) == 1
    paths_vec = node_paths[1][2]
    @assert length(paths_vec) == 1
    Ks = paths_vec[1].K
    @assert length(Ks) == 1
    return Ks[1]
end

function SOI.forward_pass!(sp::SOI.AbstractStochasticProgram, algo::Algorithm, to::TimerOutput, result::SOI.Result, verbose)
    master = SOI.get(sp, SOI.MasterNode())
    NodeT = typeof(master)
    @timeit to "solve" mastersol = SOI.get(sp, SOI.Solution(), master)
    infeasibility_detected = SOI.getstatus(mastersol) == :Infeasible

    num_stages = SOI.get(sp, SOI.NumberOfStages())
    PathT = SDDPPath{SOI.get(sp, SOI.TransitionType()), typeof(mastersol)}
    TT = Tuple{NodeT, Vector{PathT}}
    pathsd = Vector{Vector{TT}}(undef, num_stages)
    if infeasibility_detected
        pathsd[1] = TT[]
    else
        pathsd[1] = TT[(master, [SDDPPath{SOI.get(sp, SOI.TransitionType())}(mastersol, [SOI.getnodeobjectivevalue(mastersol)], [1.], [algo.K])])]
    end
    endedpaths = PathT[]

    for t in 2:num_stages
        verbose >= 3 && println("Stage $t/$num_stages")

        # Merge paths
        if algo.mergepaths
            pathsd[t-1] = mergesamepaths(pathsd[t-1], to, algo.ztol)
        end

        # Make jobs
        jobsd = childjobs(sp, pathsd[t-1], algo.pathsampler, t, num_stages, endedpaths)

        # Solve Jobs (parallelism possible here)
        infeasibility_detected |= solvejobs!(sp, jobsd, to, algo.stopatinf)

        if algo.forwardcuts
            gencuts(pathsd[t-1], sp, to, algo.ztol)
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
    result.upperbound = z_UB
    result.σ_UB = σ
    result.paths = Paths(pathsd)
    result.status = SOI.getstatus(mastersol)             # TODO status in backward pass
    result.lowerbound = SOI.getobjectivevalue(mastersol) # TODO lowerbound in backward pass
end

function SOI.backward_pass!(sp::SOI.AbstractStochasticProgram, algo::Algorithm, to::TimerOutput, result::SOI.Result, verbose)
    if algo.backwardcuts
        # We could ask a node if it has new fcuts/ocuts before resolving, that way the last stage won't be a particular case anymore
        num_stages = SOI.get(sp, SOI.NumberOfStages())
        for t in num_stages:-1:2
            # The last stages do not need to be resolved since it does not have any new cut
            if t != num_stages
                # Make jobs
                endedpaths = SDDPPath{SOI.get(sp, SOI.TransitionType())}[]
                jobsd = childjobs(sp, result.paths.paths[t-1], algo.pathsampler, t, num_stages, endedpaths) # FIXME shouldn't need pathsampler here
                # Solve Jobs (parallelism possible here)
                solvejobs!(sp, jobsd, to, algo.stopatinf)
            end

            gencuts(result.paths.paths[t-1], sp, to, algo.ztol)
            # The cut is added after so that they are added by group and duplicate detection in CutPruners works better
            applycuts(result.paths.paths[t-1], sp)
        end
    end
end
