abstract type AbstractAlgorithm end

# Solution of the full optimization, not an AbstractSolution
struct SolutionSummary
    status
    objval
    sol
    attrs
end

"""
    optimize!(sp::AbstractStochasticProgram, algo::AbstractAlgorithm, stopcrit::AbstractStoppingCriterion, verbose=0)

Run the algorithm `algo` on the stochastic program `sp` until the termination criterion `stopcrit` requires stopping with verbose level `verbose`.
"""
function optimize!(sp::AbstractStochasticProgram, algo::AbstractAlgorithm, stopcrit::AbstractStoppingCriterion=Pereira(), verbose=0)
    # Default implementation, define a specific method for algorithms for which
    # this default is not appropriate
    mastersol = nothing
    totalstats = SDDPStats()
    stats = SDDPStats()
    stats.niterations = 1

    while (mastersol === nothing || getstatus(mastersol) != :Infeasible) && !stop(stopcrit, stats, totalstats)
        itertime = @_time mastersol, stats = iterate!(sp, algo, verbose)
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
        println("Algorithm completed in $(totalstats.niterations) iterations in $(totalstats.time) s")
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
    SolutionSummary(getstatus(mastersol), getobjectivevalue(mastersol), getstatevalue(mastersol), attrs)
end

"""
    iterate!(sp::AbstractStochasticProgram, algo::AbstractAlgorithm, verbose)

Run one iteration of the algorithm `algo` on the stochastic program `sp` with verbose level `verbose`.
Return the solution of the master state and the stats.
"""
function iterate!(sp::AbstractStochasticProgram, algo::AbstractAlgorithm, verbose)
    # Default implementation, define a specific method for algorithms for which
    # this default is not appropriate
    paths, mastersol, forward_stats = forwardpass!(sp, algo, verbose)
    backward_stats = backwardpass!(sp, algo, paths, verbose)
    # Uses forward_stats in the RHS so that its values are used for upperbound, lowerbound, σ_UB and npaths
    mastersol, backward_stats + forward_stats
end

"""
    forwardpass!(sp::AbstractStochasticProgram, algo::AbstractAlgorithm, verbose)

Run the forward pass of algorithm `algo` on the stochastic program `sp` with verbose level `verbose`.
Returns the forward paths, the master state solution and the stats.
"""
function forwardpass! end

"""
    backwardpass!(sp::AbstractStochasticProgram, algo::AbstractAlgorithm, paths, verbose)

Run the backward pass of algorithm `algo` on the stochastic program `sp` for the paths `paths` with verbose level `verbose`.
Returns the stats.
"""
function backwardpass! end
