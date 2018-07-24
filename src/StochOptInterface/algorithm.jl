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
    totalstats = Stats()
    stats = Stats()
    stats.niterations = 1
    totalto = TimerOutput()
    to = totalto

    while (mastersol === nothing || getstatus(mastersol) != :Infeasible) && !stop(stopcrit, to, stats, totalstats)
        @timeit totalto "iteration $(totalstats.niterations+1)" mastersol, stats = iterate!(sp, algo, totalto, verbose)
        to = totalto["iteration $(totalstats.niterations+1)"]
        fto = TimerOutputs.flatten(totalto)
        totalstats += stats
        if verbose >= 2
            println("Stats for this iteration:")
            println(to)
            println("Status: $(getstatus(mastersol))")
            println("Solution value: $(getstatevalue(mastersol))")
            println(stats)
            println("Total stats:")
            println(fto)
            println(totalstats)
        end
    end

    fto = TimerOutputs.flatten(totalto)

    if verbose >= 1
        println(fto)
        println("Status: $(getstatus(mastersol))")
        println("Objective value: $(getobjectivevalue(mastersol))")
        println(totalstats)
    end

    attrs = Dict()
    attrs[:stats] = totalstats
    attrs[:niter] = totalstats.niterations
    attrs[:nfcuts] = nfcuts(fto)
    attrs[:nocuts] = nocuts(fto)
    SolutionSummary(getstatus(mastersol), getobjectivevalue(mastersol), getstatevalue(mastersol), attrs)
end

"""
    iterate!(sp::AbstractStochasticProgram, algo::AbstractAlgorithm, to::TimerOutput, verbose)

Run one iteration of the algorithm `algo` on the stochastic program `sp` with verbose level `verbose`.
Return the solution of the master state and the stats.
"""
function iterate!(sp::AbstractStochasticProgram, algo::AbstractAlgorithm, to::TimerOutput, verbose)
    # Default implementation, define a specific method for algorithms for which
    # this default is not appropriate
    paths, mastersol, stats = forwardpass!(sp, algo, to, verbose)
    backwardpass!(sp, algo, paths, to, verbose)
    # Uses forward_stats in the RHS so that its values are used for upperbound, lowerbound, σ_UB and npaths
    return mastersol, stats
end

"""
    forwardpass!(sp::AbstractStochasticProgram, algo::AbstractAlgorithm, to::TimerOutput, verbose)

Run the forward pass of algorithm `algo` on the stochastic program `sp` with verbose level `verbose`.
Returns the forward paths, the master state solution and the stats.
"""
function forwardpass! end

"""
    backwardpass!(sp::AbstractStochasticProgram, algo::AbstractAlgorithm, paths, to::TimerOutput, verbose)

Run the backward pass of algorithm `algo` on the stochastic program `sp` for the paths `paths` with verbose level `verbose`.
"""
function backwardpass! end
