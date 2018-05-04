export AbstractCutGenerator, AbstractOptimalityCutGenerator, NoOptimalityCutGenerator, MultiCutGenerator, AvgCutGenerator
export nθ, solveallchildren, needallchildsol
abstract type AbstractCutGenerator end
abstract type AbstractOptimalityCutGenerator <: AbstractCutGenerator end

struct FeasibilityCutGenerator <: AbstractCutGenerator
end
function gencut(::FeasibilityCutGenerator, sp, parent, path, stats, ztol)
    for tr in out_transitions(sp, parent)
        child = target(sp, tr)
        if haskey(path.childsols, child)
            childsol = path.childsols[child]
            if getstatus(childsol) == :Infeasible
                stats.nfcuts += 1
                stats.fcutstime += @_time add_feasibility_cut!(sp, child, feasibility_cut(childsol)..., parent)
            end
        end
    end
end
function applycut(::FeasibilityCutGenerator, sp, node)
    for tr in out_transitions(sp, node)
        apply_feasibility_cuts!(sp, target(sp, tr))
    end
end

# No Optimality Cut Generator
struct NoOptimalityCutGenerator <: AbstractOptimalityCutGenerator
end
nθ(::NoOptimalityCutGenerator, proba) = 0
needallchildsol(::NoOptimalityCutGenerator) = false
function gencut(::NoOptimalityCutGenerator, sp, parent, path, stats, ztol)
end
function applycut(::NoOptimalityCutGenerator, sp, node)
    apply_optimality_cuts!(sp, node)
end

# Multi Cut Generator
struct MultiCutGenerator <: AbstractOptimalityCutGenerator
end
nθ(::MultiCutGenerator, proba) = length(proba)
needallchildsol(::MultiCutGenerator) = false
function gencut(::MultiCutGenerator, sp, parent, path, stats, ztol)
    for (child, sol) in path.childsols
        status = getstatus(sol)
        if status != :Unbounded
            a, β = optimality_cut(sol)
            @assert status == :Optimal
            edge = Edge(parent, child)
            if haskey(sp.childT, edge)
                aT = sp.childT[edge]' * a
            else
                aT = a
            end
            x = getstatevalue(path.sol)
            θ = getθvalue(sp, parent, child, path.sol)
            if getstatus(path.sol) == :Unbounded || _lt(θ, β - dot(aT, x), ztol)
                stats.ocutstime += @_time add_optimality_cut_for_parent!(sp, child, a, β, parent)
                stats.nocuts += 1
            end
        end
    end
end
function applycut(::MultiCutGenerator, sp, node)
    for tr in out_transitions(sp, node)
        apply_optimality_cuts_for_parent!(sp, target(sp, tr))
    end
end

# Average Cut Generator
struct AvgCutGenerator <: AbstractOptimalityCutGenerator
end
nθ(::AvgCutGenerator, proba) = 1
needallchildsol(::AvgCutGenerator) = true
function gencut(::AvgCutGenerator, sp, parent, path, stats, ztol)
    (!path.childs_feasible || !path.childs_bounded) && return
    avga = zeros(statedim(sp, parent))
    avgβ = 0.
    for (child, sol) in path.childsols
        status = getstatus(sol)
        @assert status != :Unbounded
        a, β = optimality_cut(sol)
        @assert status == :Optimal
        edge = Edge(parent, child)
        if haskey(sp.childT, edge)
            aT = sp.childT[edge]' * a
        else
            aT = a
        end
        proba = probability(sp, edge)
        avga += proba * aT
        avgβ += proba * β
    end
    x = getstatevalue(path.sol)
    θ = getθvalue(sp, parent, path.sol)
    if getstatus(path.sol) == :Unbounded || _lt(θ, avgβ - dot(avga, x), ztol)
        stats.ocutstime += @_time add_optimality_cut!(sp, parent, avga, avgβ, parent)
        stats.nocuts += 1
    end
end
function applycut(::AvgCutGenerator, sp, node)
    apply_optimality_cuts!(sp, node)
end
