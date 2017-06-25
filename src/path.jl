function meanstdpaths(z::Vector{Float64}, proba::Vector{Float64}, npaths::Vector{Int}, Ktot)
    if Ktot != -1
        @assert sum(npaths) == Ktot
        proba = npaths / Ktot
    end
    μ = dot(proba, z)
    σ = sqrt(dot(proba, (z - μ).^2))
    μ, σ
end

type SDDPPath
    sol::NLDSSolution
    z::Vector{Float64}
    proba::Vector{Float64}
    K::Vector{Int}
    childs_feasible::Bool
    childs_bounded::Bool
    childsols::Vector{Nullable{NLDSSolution}}

    function SDDPPath(sol, z, proba, K, nchilds)
        childsols = Nullable{NLDSSolution}[nothing for i in 1:nchilds]
        new(sol, z, proba, K, true, true, childsols)
    end
end

function meanstdpaths(paths::Vector{SDDPPath}, Ktot)
    z = reduce(append!, Float64[], Vector{Float64}[x.z for x in paths])
    proba = reduce(append!, Float64[], Vector{Float64}[x.proba for x in paths])
    npaths = reduce(append!, Int[], Vector{Int}[x.K for x in paths])
    meanstdpaths(z, proba, npaths, Ktot)
end

function canmerge(p::SDDPPath, q::SDDPPath, ztol)
    _isapprox(p.sol.x, q.sol.x, ztol)
end

function merge!(p::SDDPPath, q::SDDPPath)
    @assert p.childs_feasible == q.childs_feasible
    append!(p.z, q.z)
    append!(p.proba, q.proba)
    append!(p.K, q.K)
end

type SDDPJob
    sol::Nullable{NLDSSolution}
    proba::Vector{Float64}
    K::Vector{Int}
    parentnode::SDDPNode
    parent::SDDPPath
    i::Int

    function SDDPJob(proba::Vector{Float64}, K::Vector{Int}, parentnode::SDDPNode, parent::SDDPPath, i::Int)
        new(nothing, proba, K, parentnode, parent, i::Int)
    end
end


function Base.isapprox(p::SDDPPath, q::SDDPPath)
    Base.isapprox(p.sol.x, q.sol.x)
end

function addjob!{S}(jobsd::Dict{SDDPNode{S}, Vector{SDDPJob}}, node::SDDPNode{S}, job::SDDPJob)
    if node in keys(jobsd)
        push!(jobsd[node], job)
    else
        jobsd[node] = [job]
    end
end

function mergesamepaths{StateT}(pathsd::Vector{Tuple{StateT, Vector{SDDPPath}}}, stats, ztol)
    before = sum([sum([sum(path.K) for path in paths]) for (state, paths) in pathsd])
    newpathsd = Tuple{StateT, Vector{SDDPPath}}[]
    stats.mergetime += @_time for (state, paths) in pathsd
        keep = ones(Bool, length(paths))
        merged = false
        for i in 1:length(paths)
            for j in 1:(i-1)
                if keep[j] && canmerge(paths[i], paths[j], ztol)
                    #println("Merging path since ||x_i - x_j||_∞ = $(norm(paths[j].sol.x - paths[i].sol.x, Inf))")
                    merge!(paths[i], paths[j])
                    keep[j] = false
                    merged = true
                    stats.nmerged += 1
                    break
                end
            end
        end
        push!(newpathsd, (state, paths[keep]))
    end
    after = sum([sum([sum(path.K) for path in paths]) for (state, paths) in newpathsd])
    @assert before == after
    newpathsd
end
