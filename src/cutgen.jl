export AbstractCutGenerator, NoOptimalityCutGenerator, MultiCutGenerator, AvgCutGenerator
export nθ, solveallchildren, needallchildsol
@compat abstract type AbstractCutGenerator end

# No Optimality Cut Generator
immutable NoOptimalityCutGenerator <: AbstractCutGenerator
end
nθ(::NoOptimalityCutGenerator, proba) = 0
needallchildsol(::NoOptimalityCutGenerator) = false
function gencut(::NoOptimalityCutGenerator, parent, path, stats, ztol)
end
function applycut(::NoOptimalityCutGenerator, node, g)
    applyoptimalitycut!(node)
end

# Multi Cut Generator
immutable MultiCutGenerator <: AbstractCutGenerator
end
nθ(::MultiCutGenerator, proba) = length(proba)
needallchildsol(::MultiCutGenerator) = false
function gencut(::MultiCutGenerator, parent, path, stats, ztol)
    for i in 1:length(parent.children)
        if !isnull(path.childsols[i]) && get(path.childsols[i]).status != :Unbounded
            childsol = get(path.childsols[i])
            a = childsol.πT
            β = childsol.πh + childsol.σd
            @assert childsol.status == :Optimal
            if isnull(parent.childT)
                aT = a
            else
                aT = get(parent.childT)[i]' * a
            end
            β += childsol.ρe
            if path.sol.status == :Unbounded || _lt(path.sol.θ[i], β - dot(aT, path.sol.x), ztol)
                stats.ocutstime += @_time pushoptimalitycutforparent!(parent.children[i], a, β, parent)
                stats.nocuts += 1
            end
        end
    end
end
function applycut(::MultiCutGenerator, node, g)
    for child in children(g, node)
        applyoptimalitycutforparent!(child)
    end
end

# Average Cut Generator
immutable AvgCutGenerator <: AbstractCutGenerator
end
nθ(::AvgCutGenerator, proba) = 1
needallchildsol(::AvgCutGenerator) = true
function gencut(::AvgCutGenerator, parent, path, stats, ztol)
    (!path.childs_feasible || !path.childs_bounded) && return
    avga = zeros(parent.nlds.nx)
    avgβ = 0
    for i in 1:length(parent.children)
        @assert !isnull(path.childsols[i])
        @assert get(path.childsols[i]).status != :Unbounded
        childsol = get(path.childsols[i])
        a = childsol.πT
        β = childsol.πh + childsol.σd
        @assert childsol.status == :Optimal
        if isnull(parent.childT)
            aT = a
        else
            aT = get(parent.childT)[i]' * a
        end
        β += childsol.ρe
        avga += parent.proba[i] * aT
        avgβ += parent.proba[i] * β
    end
    if path.sol.status == :Unbounded || _lt(path.sol.θ[1], avgβ - dot(avga, path.sol.x), ztol)
        stats.ocutstime += @_time pushoptimalitycut!(parent, avga, avgβ, parent)
        stats.nocuts += 1
    end
end
function applycut(::AvgCutGenerator, node, g)
    applyoptimalitycut!(node)
end
