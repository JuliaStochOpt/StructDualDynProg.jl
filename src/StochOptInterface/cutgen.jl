export AbstractCutGenerator, AbstractOptimalityCutGenerator, NoOptimalityCutGenerator, MultiCutGenerator, AvgCutGenerator
export nθ, solveallchildren, needallchildsol

abstract type AbstractCutGenerator end
abstract type AbstractCut end
abstract type AbstractOptimalityCutGenerator <: AbstractCutGenerator end
abstract type AbstractOptimalityCut <: AbstractCut end

"""
    add_cut!(sp::AbstractStochasticProgram, tr::AbstractTransition, cut::AbstractCut)

At the cut `cut` to the transition `tr`.
"""
function add_cut! end

function add_feasibility_cut! end
function add_optimality_cut! end
function add_optimality_cut_for_parent! end
function apply_feasibility_cuts! end
function apply_optimality_cuts! end
function apply_optimality_cuts_for_parent! end

struct FeasibilityCutGenerator <: AbstractCutGenerator end
#struct FeasibilityCut{T, VT::AbstractVector{T}} <: AbstractCut
#    a::VT
#    β::T
#end
function gencut(::FeasibilityCutGenerator, sp::AbstractStochasticProgram, parent, pool::AbstractSolutionPool, stats, ztol)
    for tr in get(sp, OutTransitions(), parent)
        if hassolution(pool, tr)
            trsol = getsolution(pool, tr)
            if getstatus(trsol) == :Infeasible
                stats.nfcuts += 1
                stats.fcutstime += @_time add_feasibility_cut!(sp, get(sp, Target(), tr), feasibility_cut(trsol)..., parent)
            end
        end
    end
end
function applycut(::FeasibilityCutGenerator, sp::AbstractStochasticProgram, node)
    for tr in get(sp, OutTransitions(), node)
        apply_feasibility_cuts!(sp, get(sp, Target(), tr))
    end
end

# No Optimality Cut Generator
struct NoOptimalityCutGenerator <: AbstractOptimalityCutGenerator
end
nθ(::NoOptimalityCutGenerator, proba) = 0
needallchildsol(::NoOptimalityCutGenerator) = false
function gencut(::NoOptimalityCutGenerator, sp::AbstractStochasticProgram, parent, pool::AbstractSolutionPool, stats, ztol)
end
function applycut(::NoOptimalityCutGenerator, sp::AbstractStochasticProgram, node)
    apply_optimality_cuts!(sp, node)
end

# Multi Cut Generator
struct MultiCutGenerator <: AbstractOptimalityCutGenerator
end
nθ(::MultiCutGenerator, proba) = length(proba)
needallchildsol(::MultiCutGenerator) = false
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
                    stats.ocutstime += @_time add_optimality_cut_for_parent!(sp, get(sp, Target(), tr), a, β, parent)
                    stats.nocuts += 1
                end
            end
        end
    end
end
function applycut(::MultiCutGenerator, sp::AbstractStochasticProgram, node)
    for tr in get(sp, OutTransitions(), node)
        apply_optimality_cuts_for_parent!(sp, get(sp, Target(), tr))
    end
end

# Average Cut Generator
struct AvgCutGenerator <: AbstractOptimalityCutGenerator
end
nθ(::AvgCutGenerator, proba) = 1
needallchildsol(::AvgCutGenerator) = true
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
        stats.ocutstime += @_time add_optimality_cut!(sp, parent, avga, avgβ, parent)
        stats.nocuts += 1
    end
end
function applycut(::AvgCutGenerator, sp::AbstractStochasticProgram, node)
    apply_optimality_cuts!(sp, node)
end
