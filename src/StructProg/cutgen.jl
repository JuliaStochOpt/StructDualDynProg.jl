export MultiCutGenerator, AvgCutGenerator

abstract type AbstractCutGenerator end
abstract type AbstractOptimalityCutGenerator <: AbstractCutGenerator end

struct FeasibilityCutGenerator <: AbstractCutGenerator end
function gencut(::FeasibilityCutGenerator, sp::SOI.AbstractStochasticProgram, parent, pool::SOI.AbstractSolutionPool, to::TimerOutput, ztol)
    for tr in SOI.get(sp, SOI.OutTransitions(), parent)
        if SOI.hassolution(pool, tr)
            trsol = SOI.getsolution(pool, tr)
            if SOI.getstatus(trsol) == :Infeasible
                @timeit to SOI.FCUTS_KEY SOI.addcut!(sp, tr, FeasibilityCut(SOI.feasibility_cut(trsol)...))
            end
        end
    end
end
function applycut(::FeasibilityCutGenerator, sp::SOI.AbstractStochasticProgram, node)
    for tr in SOI.get(sp, SOI.OutTransitions(), node)
        SOI.applycuts!(sp, tr, FeasibilityCut)
    end
end

# No Optimality Cut Generator
struct NoOptimalityCutGenerator <: AbstractOptimalityCutGenerator
end
nθ(::NoOptimalityCutGenerator, proba) = 0
needallsolutions(::NoOptimalityCutGenerator) = false
function gencut(::NoOptimalityCutGenerator, sp::SOI.AbstractStochasticProgram, parent, pool::SOI.AbstractSolutionPool, to::TimerOutput, ztol)
end
function applycut(::NoOptimalityCutGenerator, sp::SOI.AbstractStochasticProgram, node)
end

# Multi Cut Generator
struct MultiCutGenerator <: AbstractOptimalityCutGenerator
end
nθ(::MultiCutGenerator, proba) = length(proba)
needallsolutions(::MultiCutGenerator) = false
function gencut(::MultiCutGenerator, sp::SOI.AbstractStochasticProgram, parent, pool::SOI.AbstractSolutionPool, to::TimerOutput, ztol)
    for tr in SOI.get(sp, SOI.OutTransitions(), parent)
        if SOI.hassolution(pool, tr)
            trsol = SOI.getsolution(pool, tr)
            status = SOI.getstatus(trsol)
            if status != :Unbounded
                a, β = SOI.optimality_cut(trsol)
                @assert status == :Optimal
                if tr.childT !== nothing
                    aT = tr.childT' * a
                else
                    aT = a
                end
                sol = SOI.getsolution(pool)
                x = SOI.getnodevalue(sol)
                θ = SOI.getbellmanvalue(sp, tr, sol)
                if SOI.getstatus(sol) == :Unbounded || _lt(θ, β - dot(aT, x), ztol)
                    @timeit to SOI.OCUTS_KEY SOI.addcut!(sp, tr, MultiOptimalityCut(a, β))
                end
            end
        end
    end
end
function applycut(::MultiCutGenerator, sp::SOI.AbstractStochasticProgram, node)
    for tr in SOI.get(sp, SOI.OutTransitions(), node)
        SOI.applycuts!(sp, tr, MultiOptimalityCut)
    end
end

# Average Cut Generator
struct AvgCutGenerator <: AbstractOptimalityCutGenerator
end
nθ(::AvgCutGenerator, proba) = 1
needallsolutions(::AvgCutGenerator) = true
function gencut(::AvgCutGenerator, sp::SOI.AbstractStochasticProgram, parent, pool::SOI.AbstractSolutionPool, to::TimerOutput, ztol)
    # We need all transitions to be solved, feasible and bounded to generate an averaged cut
    (!SOI.allfeasible(pool) || !SOI.allbounded(pool)) && return
    avga = zeros(SOI.get(sp, SOI.Dimension(), parent))
    avgβ = 0.
    for tr in SOI.get(sp, SOI.OutTransitions(), parent)
        SOI.hassolution(pool, tr) || error("Average Cut Generator needs a solution for each transition")
        trsol = SOI.getsolution(pool, tr)
        status = SOI.getstatus(trsol)
        @assert status != :Unbounded
        a, β = SOI.optimality_cut(trsol)
        @assert status == :Optimal
        if tr.childT !== nothing
            aT = tr.childT' * a
        else
            aT = a
        end
        proba = SOI.get(sp, SOI.Probability(), tr)
        avga += proba * aT
        avgβ += proba * β
    end
    sol = SOI.getsolution(pool)
    x = SOI.getnodevalue(sol)
    θ = SOI.getbellmanvalue(sp, parent, sol)
    if SOI.getstatus(sol) == :Unbounded || _lt(θ, avgβ - dot(avga, x), ztol)
        @timeit to SOI.OCUTS_KEY SOI.addcut!(sp, parent, AveragedOptimalityCut(avga, avgβ))
    end
end
function applycut(::AvgCutGenerator, sp::SOI.AbstractStochasticProgram, node)
    SOI.applycuts!(sp, node, AveragedOptimalityCut)
end
