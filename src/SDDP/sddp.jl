export SDDP

struct SDDPSolution
    status
    objval
    sol
    attrs
end

function solvejobs!(sp, jobsd, stats, stopatinf)
    infeasibility_detected = false
    for (node, jobs) in jobsd
        for job in jobs
            if !stopatinf || job.parent.childs_feasible
                solvejob!(sp, job, node, stats)
                if !job.parent.childs_feasible
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
            if path.childs_feasible
                gencut(cutgenerator(sp, parent), sp, parent, path, stats, ztol)
            else
                gencut(FeasibilityCutGenerator(), sp, parent, path, stats, ztol)
            end
        end
    end
end

function applycuts(pathsd, sp)
    for (node, _) in pathsd
        applycut(FeasibilityCutGenerator(), sp, node)
        applycut(cutgenerator(sp, node), sp, node)
    end
end

function forwardpass!(sp::AbstractStochasticProgram, Ktot::Int, num_stages, verbose, pathsampler; ztol=1e-6, stopatinf=false, mergepaths=true, forwardcuts=false)
    stats = SDDPStats()

    master = getmaster(sp)
    NodeT = typeof(master)
    stats.solvertime += @_time mastersol = solve!(sp, master)
    stats.nsolved += 1
    stats.niterations += 1
    infeasibility_detected = getstatus(mastersol) == :Infeasible

    PathT = SDDPPath{transitiontype(sp), typeof(mastersol)}
    TT = Tuple{NodeT, Vector{PathT}}
    pathsd = Vector{Vector{TT}}(num_stages)
    if infeasibility_detected
        pathsd[1] = TT[]
    else
        pathsd[1] = TT[(master, [SDDPPath{transitiontype(sp)}(mastersol, [getstateobjectivevalue(mastersol)], [1.], [Ktot], outdegree(sp, master))])]
    end
    endedpaths = PathT[]

    for t in 2:num_stages
        verbose >= 3 && println("Stage $t/$num_stages")

        # Merge paths
        if mergepaths
            pathsd[t-1] = mergesamepaths(pathsd[t-1], stats, ztol)
        end

        # Make jobs
        jobsd = childjobs(sp, pathsd[t-1], pathsampler, t, num_stages, endedpaths)

        # Solve Jobs (parallelism possible here)
        infeasibility_detected |= solvejobs!(sp, jobsd, stats, stopatinf)

        if forwardcuts
            gencuts(pathsd[t-1], sp, stats, ztol)
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
        z_UB, σ = meanstdpaths(endedpaths, Ktot)
    end

    # update stats
    stats.upperbound = z_UB
    stats.σ_UB = σ
    stats.npaths = Ktot
    stats.lowerbound = getobjectivevalue(mastersol)

    pathsd, mastersol, stats
end

function backwardpass!(sp::AbstractStochasticProgram, num_stages, pathsd, pathsampler, stats; ztol=1e-6, stopatinf=false)
    # We could ask a node if it has new fcuts/ocuts before resolving, that way the last stage won't be a particular case anymore
    for t in num_stages:-1:2
        # The last stages do not need to be resolved since it does not have any new cut
        if t != num_stages
            # Make jobs
            endedpaths = SDDPPath{transitiontype(sp)}[]
            jobsd = childjobs(sp, pathsd[t-1], pathsampler, t, num_stages, endedpaths) # FIXME shouldn't need pathsampler here
            # Solve Jobs (parallelism possible here)
            solvejobs!(sp, jobsd, stats, stopatinf)
        end

        gencuts(pathsd[t-1], sp, stats, ztol)
        # The cut is added after so that they are added by group and duplicate detection in CutPruners works better
        applycuts(pathsd[t-1], sp)
    end
end

"""
$(SIGNATURES)

Runs one iteration of the SDDP algorithm on the stochastic program given by `sp`.
A total of `Ktot` paths will be explored up to `num_stages` stages.
The paths will be selected according to `pathsampler` and equivalent paths might be merged if their difference is smaller than `ztol` and `mergepaths` is true.
The parameter `ztol` is also used to check whether a new cut is useful.
When a scenario is infeasible and `stopatinf` is true then no other scenario with the same ancestor is run. Note that since the order in which the different scenarios is run is not deterministic, this might introduce nondeterminism even if the sampling is deterministic.
"""
function iteration!(sp::AbstractStochasticProgram, Ktot::Int, num_stages, verbose, pathsampler; ztol=1e-6, stopatinf=false, mergepaths=true, forwardcuts=false, backwardcuts=true)

    pathsd, mastersol, stats = forwardpass!(sp, Ktot, num_stages, verbose, pathsampler; ztol=ztol, stopatinf=stopatinf, mergepaths=mergepaths, forwardcuts=forwardcuts)

    if backwardcuts
        backwardpass!(sp, num_stages, pathsd, pathsampler, stats; ztol=ztol)
    end

    mastersol, stats
end

"""
$(SIGNATURES)

Runs the SDDP algorithms on the stochastic program given by `sp`.
The algorithm will do iterations until `stopcrit` decides to stop or when the root node is infeasible.
In each iterations, `K` paths will be explored up to `num_stages` stages.
The paths will be selected according to `pathsampler` and equivalent paths might be merged if their difference is smaller than `ztol` and `mergepaths` is true.
The parameter `ztol` is also used to check whether a new cut is useful.
When a scenario is infeasible and `stopatinf` is true then no other scenario with the same ancestor is run. Note that since the order in which the different scenarios is run is not deterministic, this might introduce nondeterminism even if the sampling is deterministic.
By default, the cuts are added backward. However, if `forwardcuts` is set to `true` and `backwardcuts` is set to `false` the cuts are added forward.
"""
function SDDP(sp::AbstractStochasticProgram, num_stages; K::Int=25, stopcrit::AbstractStoppingCriterion=Pereira(), verbose=0, pathsampler::AbstractPathSampler=ProbaPathSampler(true), ztol=1e-6, stopatinf=false, mergepaths=true, forwardcuts=false, backwardcuts=true)
    mastersol = nothing
    totalstats = SDDPStats()
    stats = SDDPStats()
    stats.niterations = 1

    while (mastersol === nothing || getstatus(mastersol) != :Infeasible) && !stop(stopcrit, stats, totalstats)
        itertime = @_time mastersol, stats = iteration!(sp, K, num_stages, verbose, pathsampler, ztol=1e-6, stopatinf=stopatinf, mergepaths=mergepaths, forwardcuts=forwardcuts, backwardcuts=backwardcuts)
        stats.time = itertime

        totalstats += stats
        if verbose >= 2
            println("Iteration $(totalstats.niterations) completed in $itertime s (Total time is $(totalstats.time))")
            println("Status: $(getstatus(mastersol))")
            println("Upper Bound: $(totalstats.upperbound)")
            println("Lower Bound: $(totalstats.lowerbound)")
            #println(" Solution value: $(getstatevalue(mastersol))")
            println("Stats for this iteration:")
            println(stats)
            println("Total stats:")
            println(totalstats)
        end
    end

    if verbose >= 1
        println("SDDP completed in $(totalstats.niterations) iterations in $(totalstats.time) s")
        println("Status: $(getstatus(mastersol))")
        #println("Objective value: $(getobjectivevalue(mastersol))")
        println("Upper Bound: $(totalstats.upperbound)")
        println("Lower Bound: $(totalstats.lowerbound)")
        #println(" Solution value: $(getstatevalue(mastersol))")
        println("Total stats:")
        println(totalstats)
    end

    attrs = Dict()
    attrs[:stats] = totalstats
    attrs[:niter] = totalstats.niterations
    attrs[:nfcuts] = totalstats.nfcuts
    attrs[:nocuts] = totalstats.nocuts
    SDDPSolution(getstatus(mastersol), getobjectivevalue(mastersol), getstatevalue(mastersol), attrs)
end
