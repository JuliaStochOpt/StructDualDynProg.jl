function meanstdpaths(z::Vector{Float64}, proba::Vector{Float64}, npaths::Vector{Int}, Ktot)
    if Ktot != -1
        @assert sum(npaths) == Ktot
        proba = npaths / Ktot
    end
    μ = dot(proba, z)
    σ = sqrt(dot(proba, (z - μ).^2))
    μ, σ
end

mutable struct SDDPPath
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

type SDDPJob{NodeT}
    sol::Nullable{NLDSSolution}
    proba::Vector{Float64}
    K::Vector{Int}
    parentnode::NodeT
    parent::SDDPPath
    i::Int

    function SDDPJob{NodeT}(proba::Vector{Float64}, K::Vector{Int}, parentnode::NodeT, parent::SDDPPath, i::Int) where {NodeT}
        new(nothing, proba, K, parentnode, parent, i::Int)
    end
end

function Base.isapprox(p::SDDPPath, q::SDDPPath)
    Base.isapprox(p.sol.x, q.sol.x)
end

function addjob!{NodeT}(jobsd::Dict{NodeT, Vector{SDDPJob{NodeT}}}, node::NodeT, job::SDDPJob{NodeT})
    if node in keys(jobsd)
        push!(jobsd[node], job)
    else
        jobsd[node] = [job]
    end
end

function mergesamepaths{NodeT}(pathsd::Vector{Tuple{NodeT, Vector{SDDPPath}}}, stats, ztol)
    before = sum([sum([sum(path.K) for path in paths]) for (node, paths) in pathsd])
    newpathsd = Tuple{NodeT, Vector{SDDPPath}}[]
    stats.mergetime += @_time for (node, paths) in pathsd
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
        push!(newpathsd, (node, paths[keep]))
    end
    after = sum([sum([sum(path.K) for path in paths]) for (node, paths) in newpathsd])
    @assert before == after
    newpathsd
end

function childjobs{NodeT}(g::AbstractSDDPGraph, pathsd::Vector{Tuple{NodeT, Vector{SDDPPath}}}, pathsampler, t, num_stages)
    jobsd = Dict{NodeT, Vector{SDDPJob{NodeT}}}()
    for (node, paths) in pathsd
        if haschildren(g, node)
            for path in paths
                # Adding Jobs
                npaths = samplepaths(pathsampler, g, node, path.K, t, num_stages)
                childocuts = Array{Any}(nchildren(g, node))
                for i in 1:nchildren(g, node)
                    if sum(npaths[i]) != 0 || needallchildsol(cutgen(g, node)) # || t == 2
                        addjob!(jobsd, getchild(g, node, i), SDDPJob{NodeT}(path.proba * getproba(g, node, i), npaths[i], node, path, i))
                    end
                end
            end
        else
            append!(endedpaths, paths)
        end
    end
    jobsd
end

function jobstopaths{NodeT}(jobsd::Dict{NodeT, Vector{SDDPJob{NodeT}}}, g::AbstractSDDPGraph)
    pathsd = Tuple{NodeT, Vector{SDDPPath}}[]
    for (node, jobs) in jobsd
        K = [find(job.K .!= 0) for job in jobs]
        keep = find(Bool[jobs[i].parent.childs_feasible && !isempty(K[i]) for i in 1:length(jobs)])
        if !isempty(keep)
            paths = SDDPPath[SDDPPath(get(jobs[i].sol), jobs[i].parent.z[K[i]]+get(jobs[i].sol).objvalx, jobs[i].proba[K[i]], jobs[i].K[K[i]], nchildren(g, node)) for i in keep]
            push!(pathsd, (node, paths))
        end
    end
    pathsd
end

function solvejob!(job::SDDPJob, node, stats)
    stats.setxtime += @_time setchildx(job.parentnode, job.i, job.parent.sol)
    stats.nsetx += 1
    stats.solvertime += @_time job.sol = loadAndSolve(node)
    job.parent.childsols[job.i] = job.sol
    stats.nsolved += 1
    if get(job.sol).status == :Infeasible
        job.parent.childs_feasible = false
    elseif get(job.sol).status == :Unbounded
        job.parent.childs_bounded = false
    end
end
