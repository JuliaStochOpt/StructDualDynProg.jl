export AbstractCutGenerator, AbstractOptimalityCutGenerator, NoOptimalityCutGenerator, MultiCutGenerator, AvgCutGenerator

abstract type AbstractCutGenerator end
abstract type AbstractOptimalityCutGenerator <: AbstractCutGenerator end

struct FeasibilityCutGenerator <: AbstractCutGenerator end
function gencut(::FeasibilityCutGenerator, sp::AbstractStochasticProgram, parent, pool::AbstractSolutionPool, stats, ztol)
    for tr in get(sp, OutTransitions(), parent)
        if hassolution(pool, tr)
            trsol = getsolution(pool, tr)
            if getstatus(trsol) == :Infeasible
                stats.nfcuts += 1
                stats.fcutstime += @_time addcut!(sp, tr, FeasibilityCut(feasibility_cut(trsol)...))
            end
        end
    end
end
function applycut(::FeasibilityCutGenerator, sp::AbstractStochasticProgram, node)
    for tr in get(sp, OutTransitions(), node)
        applycuts!(sp, tr, FeasibilityCut)
    end
end

# No Optimality Cut Generator
struct NoOptimalityCutGenerator <: AbstractOptimalityCutGenerator
end
nθ(::NoOptimalityCutGenerator, proba) = 0
needallsolutions(::NoOptimalityCutGenerator) = false
function gencut(::NoOptimalityCutGenerator, sp::AbstractStochasticProgram, parent, pool::AbstractSolutionPool, stats, ztol)
end
function applycut(::NoOptimalityCutGenerator, sp::AbstractStochasticProgram, node)
end

# Multi Cut Generator
struct MultiCutGenerator <: AbstractOptimalityCutGenerator
end
nθ(::MultiCutGenerator, proba) = length(proba)
needallsolutions(::MultiCutGenerator) = false
function gencut(::MultiCutGenerator, sp::AbstractStochasticProgram, parent, pool::AbstractSolutionPool, stats, ztol)
    for tr in get(sp, OutTransitions(), parent)
        if hassolution(pool, tr)
            trsol = getsolution(pool, tr)
            status = getstatus(trsol)
            if status != :Unbounded
                a, β = optimality_cut(trsol)
                @assert status == :Optimal
                if !isnull(tr.childT)
                    aT = get(tr.childT)' * a
                else
                    aT = a
                end
                sol = getsolution(pool)
                x = getstatevalue(sol)
                θ = getθvalue(sp, tr, sol)
                if getstatus(sol) == :Unbounded || _lt(θ, β - dot(aT, x), ztol)
                    stats.ocutstime += @_time addcut!(sp, tr, MultiOptimalityCut(a, β))
                    stats.nocuts += 1
                end
            end
        end
    end
end
function applycut(::MultiCutGenerator, sp::AbstractStochasticProgram, node)
    for tr in get(sp, OutTransitions(), node)
        applycuts!(sp, tr, MultiOptimalityCut)
    end
end

# Average Cut Generator
struct AvgCutGenerator <: AbstractOptimalityCutGenerator
end
nθ(::AvgCutGenerator, proba) = 1
needallsolutions(::AvgCutGenerator) = true
function gencut(::AvgCutGenerator, sp::AbstractStochasticProgram, parent, pool::AbstractSolutionPool, stats, ztol)
    # We need all transitions to be solved, feasible and bounded to generate an averaged cut
    (!allfeasible(pool) || !allbounded(pool)) && return
    avga = zeros(get(sp, Dimension(), parent))
    avgβ = 0.
    for tr in get(sp, OutTransitions(), parent)
        hassolution(pool, tr) || error("Average Cut Generator needs a solution for each transition")
        trsol = getsolution(pool, tr)
        status = getstatus(trsol)
        @assert status != :Unbounded
        a, β = optimality_cut(trsol)
        @assert status == :Optimal
        if !isnull(tr.childT)
            aT = get(tr.childT)' * a
        else
            aT = a
        end
        proba = get(sp, Probability(), tr)
        avga += proba * aT
        avgβ += proba * β
    end
    sol = getsolution(pool)
    x = getstatevalue(sol)
    θ = getθvalue(sp, parent, sol)
    if getstatus(sol) == :Unbounded || _lt(θ, avgβ - dot(avga, x), ztol)
        stats.ocutstime += @_time addcut!(sp, parent, AveragedOptimalityCut(avga, avgβ))
        stats.nocuts += 1
    end
end
function applycut(::AvgCutGenerator, sp::AbstractStochasticProgram, node)
    applycuts!(sp, node, AveragedOptimalityCut)
end
