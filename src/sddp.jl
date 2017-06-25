export SDDP

struct SDDPSolution
    status
    objval
    sol
    attrs
end

function solvejobs!(jobsd, stats, stopatinf)
    infeasibility_detected = false
    for (node, jobs) in jobsd
        for job in jobs
            if !stopatinf || job.parent.childs_feasible
                solvejob!(job, node, stats)
                if !job.parent.childs_feasible
                    infeasibility_detected = true
                end
            end
        end
    end
    infeasibility_detected
end

function gencuts(pathsd, g, stats, ztol)
    for (parent, paths) in pathsd
        for path in paths
            if path.childs_feasible
                gencut(cutgen(g, parent), g, parent, path, stats, ztol)
            else
                gencut(FeasibilityCutGenerator(), g, parent, path, stats, ztol)
            end
        end
    end
end

function applycuts(pathsd, g)
    for (node, _) in pathsd
        applycut(FeasibilityCutGenerator(), g, node)
        applycut(cutgen(g, node), g, node)
    end
end

"""
$(SIGNATURES)

Runs one iteration of the SDDP algorithm on the lattice given by `g`.
A total of `Ktot` paths will be explored up to `num_stages` stages.
The paths will be selected according to `pathsampler` and equivalent paths might be merged if their difference is smaller than `ztol` and `mergepaths` is true.
The parameter `ztol` is also used to check whether a new cut is useful.
When a scenario is infeasible and `stopatinf` is true then no other scenario with the same ancestor is run. Note that since the order in which the different scenarios is run is not deterministic, this might introduce nondeterminism even if the sampling is deterministic.
"""
function iteration{S}(g::AbstractSDDPGraph{S}, Ktot::Int, num_stages, verbose, pathsampler; ztol=1e-6, stopatinf=false, mergepaths=true, forwardcuts=false, backwardcuts=true)
    stats = SDDPStats()

	master, initialnode = getmaster(g)
	NodeT = typeof(initialnode)
    stats.solvertime += @_time mastersol = loadAndSolve(master)
    stats.nsolved += 1
    stats.niterations += 1
    infeasibility_detected = mastersol.status == :Infeasible

    pathsd = Vector{Vector{Tuple{NodeT, Vector{SDDPPath}}}}(num_stages)
    if infeasibility_detected
        pathsd[1] = Tuple{NodeT, Vector{SDDPPath}}[]
    else
        pathsd[1] = Tuple{NodeT, Vector{SDDPPath}}[(initialnode, [SDDPPath(mastersol, [mastersol.objvalx], [1.], [Ktot], length(master.children))])]
    end
    endedpaths = SDDPPath[]

    for t in 2:num_stages
        verbose >= 3 && println("Stage $t/$num_stages")

        # Merge paths
        if mergepaths
            pathsd[t-1] = mergesamepaths(pathsd[t-1], stats, ztol)
        end

        # Make jobs
        jobsd = childjobs(g, pathsd[t-1], pathsampler, t, num_stages)

        # Solve Jobs (parallelism possible here)
        infeasibility_detected |= solvejobs!(jobsd, stats, stopatinf)

        if forwardcuts
            gencuts(pathsd[t-1], g, stats, ztol)
            # The cut is added after so that they are added by group and duplicate detection in CutPruners works better
            applycuts(pathsd[t-1], g)
        end

        # Jobs -> Paths
        pathsd[t] = jobstopaths(jobsd, g)
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

    if backwardcuts # We could ask a node if it has new fcuts/ocuts before resolving, that way the last stage won't be a particular case anymore
        for t in num_stages:-1:2
            # The last stages do not need to be resolved since it does not have any new cut
            if t != num_stages
                # Make jobs
                jobsd = childjobs(g, pathsd[t-1], pathsampler, t, num_stages)
                # Solve Jobs (parallelism possible here)
                infeasibility_detected |= solvejobs!(jobsd, stats, stopatinf)
            end

            gencuts(pathsd[t-1], g, stats, ztol)
            # The cut is added after so that they are added by group and duplicate detection in CutPruners works better
            applycuts(pathsd[t-1], g)
        end
    end

    # update stats
    stats.upperbound = z_UB
    stats.σ_UB = σ
    stats.npaths = Ktot
    stats.lowerbound = mastersol.objval

    mastersol, stats
end

"""
$(SIGNATURES)

Runs the SDDP algorithms on the lattice given by `g`.
The algorithm will do iterations until `stopcrit` decides to stop or when the root node is infeasible.
In each iterations, `K` paths will be explored up to `num_stages` stages.
The paths will be selected according to `pathsampler` and equivalent paths might be merged if their difference is smaller than `ztol` and `mergepaths` is true.
The parameter `ztol` is also used to check whether a new cut is useful.
When a scenario is infeasible and `stopatinf` is true then no other scenario with the same ancestor is run. Note that since the order in which the different scenarios is run is not deterministic, this might introduce nondeterminism even if the sampling is deterministic.
By default, the cuts are added backward. However, if `forwardcuts` is set to `true` and `backwardcuts` is set to `false` the cuts are added forward.
"""
function SDDP(g::AbstractSDDPGraph, num_stages; K::Int=25, stopcrit::AbstractStoppingCriterion=Pereira(), verbose=0, pathsampler::AbstractPathSampler=ProbaPathSampler(true), ztol=1e-6, stopatinf=false, mergepaths=true, forwardcuts=false, backwardcuts=true)
    mastersol = nothing
    totalstats = SDDPStats()
    stats = SDDPStats()
    stats.niterations = 1

    while (mastersol === nothing || mastersol.status != :Infeasible) && !stop(stopcrit, stats, totalstats)
        itertime = @_time mastersol, stats = iteration(g, K, num_stages, verbose, pathsampler, ztol=1e-6, stopatinf=stopatinf, mergepaths=mergepaths, forwardcuts=forwardcuts, backwardcuts=backwardcuts)
        stats.time = itertime

        totalstats += stats
        if verbose >= 2
            println("Iteration $(totalstats.niterations) completed in $itertime s (Total time is $(stats.time))")
            println("Status: $(mastersol.status)")
            println("Upper Bound: $(totalstats.upperbound)")
            println("Lower Bound: $(totalstats.lowerbound)")
            #println(" Solution value: $(mastersol.x)")
            println("Stats for this iteration:")
            println(stats)
            println("Total stats:")
            println(totalstats)
        end
    end

    if verbose >= 1
        println("SDDP completed in $(totalstats.niterations) iterations in $(totalstats.time) s")
        println("Status: $(mastersol.status)")
        #println("Objective value: $(mastersol.objval)")
        println("Upper Bound: $(totalstats.upperbound)")
        println("Lower Bound: $(totalstats.lowerbound)")
        #println(" Solution value: $(mastersol.x)")
        println("Total stats:")
        println(totalstats)
    end

    attrs = Dict()
    attrs[:stats] = totalstats
    attrs[:niter] = totalstats.niterations
    attrs[:nfcuts] = totalstats.nfcuts
    attrs[:nocuts] = totalstats.nocuts
    SDDPSolution(mastersol.status, mastersol.objval, mastersol.x, attrs)
end
