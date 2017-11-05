export AbstractCutGenerator, AbstractOptimalityCutGenerator, NoOptimalityCutGenerator, MultiCutGenerator, AvgCutGenerator
export nθ, solveallchildren, needallchildsol
abstract type AbstractCutGenerator end
abstract type AbstractOptimalityCutGenerator <: AbstractCutGenerator end

struct FeasibilityCutGenerator <: AbstractCutGenerator
end
function gencut(::FeasibilityCutGenerator, sp, parent, path, stats, ztol)
    for child in out_neighbors(sp, parent)
        if haskey(path.childsols, child)
            childsol = path.childsols[child]
            if childsol.status == :Infeasible
                stats.nfcuts += 1
                stats.fcutstime += @_time add_feasibility_cut!(sp, child, getfeasibilitycut(childsol)..., parent)
            end
        end
    end
end
function applycut(::FeasibilityCutGenerator, sp, node)
    for child in out_neighbors(sp, node)
        apply_feasibility_cuts!(sp, child)
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
        if sol.status != :Unbounded
            a, β = getoptimalitycut(sol)
            @assert sol.status == :Optimal
            edge = Edge(parent, child)
            if haskey(sp.childT, edge)
                aT = sp.childT[edge]' * a
            else
                aT = a
            end
            θ = path.sol.θ[edgeid(sp, Edge(parent, child))]
            if path.sol.status == :Unbounded || _lt(θ, β - dot(aT, path.sol.x), ztol)
                stats.ocutstime += @_time add_optimality_cut_for_parent!(sp, child, a, β, parent)
                stats.nocuts += 1
            end
        end
    end
end
function applycut(::MultiCutGenerator, sp, node)
    for child in out_neighbors(sp, node)
        apply_optimality_cuts_for_parent!(sp, child)
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
        @assert sol.status != :Unbounded
        a, β = getoptimalitycut(sol)
        @assert sol.status == :Optimal
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
    if path.sol.status == :Unbounded || _lt(path.sol.θ[1], avgβ - dot(avga, path.sol.x), ztol)
        stats.ocutstime += @_time add_optimality_cut!(sp, parent, avga, avgβ, parent)
        stats.nocuts += 1
    end
end
function applycut(::AvgCutGenerator, sp, node)
    apply_optimality_cuts!(sp, node)
end
