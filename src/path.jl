function meanstdpaths(z::Vector{Float64}, proba::Vector{Float64}, npaths::Vector{Int}, Ktot)
    if Ktot != -1
        @assert sum(npaths) == Ktot
        proba = npaths / Ktot
    end
    μ = dot(proba, z)
    σ = sqrt(dot(proba, (z - μ).^2))
    μ, σ
end

mutable struct SDDPPath{SolT}
    sol::SolT
    z::Vector{Float64}
    proba::Vector{Float64}
    K::Vector{Int}
    childs_feasible::Bool
    childs_bounded::Bool
    childsols::Dict{Int, SolT}

    function SDDPPath(sol::SolT, z, proba, K, nchilds) where SolT
        new{SolT}(sol, z, proba, K, true, true, Dict{Int, SolT}())
    end
end

function meanstdpaths(paths::Vector{<:SDDPPath}, Ktot)
    z = reduce(append!, Float64[], Vector{Float64}[x.z for x in paths])
    proba = reduce(append!, Float64[], Vector{Float64}[x.proba for x in paths])
    npaths = reduce(append!, Int[], Vector{Int}[x.K for x in paths])
    meanstdpaths(z, proba, npaths, Ktot)
end

function Base.isapprox(p::SDDPPath, q::SDDPPath)
    Base.isapprox(getstatevalue(p.sol), getstatevalue(q.sol))
end

function canmerge(p::SDDPPath, q::SDDPPath, ztol)
    _isapprox(getstatevalue(p.sol), getstatevalue(q.sol), ztol)
end

function merge!(p::SDDPPath, q::SDDPPath)
    @assert p.childs_feasible == q.childs_feasible
    append!(p.z, q.z)
    append!(p.proba, q.proba)
    append!(p.K, q.K)
end

mutable struct SDDPJob{SolT, NodeT}
    sol::Nullable{SolT}
    proba::Vector{Float64}
    K::Vector{Int}
    parentnode::NodeT
    parent::SDDPPath{SolT}
    child::NodeT

    function SDDPJob{NodeT}(proba::Vector{Float64}, K::Vector{Int}, parentnode::NodeT, parent::SDDPPath{SolT}, child::NodeT) where {SolT, NodeT}
        new{SolT, NodeT}(nothing, proba, K, parentnode, parent, child)
    end
end

function addjob!(jobsd::Dict{NodeT, Vector{SDDPJob{SolT, NodeT}}}, node::NodeT, job::SDDPJob{SolT, NodeT}) where {SolT, NodeT}
    if node in keys(jobsd)
        push!(jobsd[node], job)
    else
        jobsd[node] = [job]
    end
end

function mergesamepaths(pathsd::Vector{Tuple{NodeT, Vector{SDDPPath{SolT}}}}, stats, ztol) where {SolT, NodeT}
    before = sum([sum([sum(path.K) for path in paths]) for (node, paths) in pathsd])
    newpathsd = Tuple{NodeT, Vector{SDDPPath{SolT}}}[]
    stats.mergetime += @_time for (node, paths) in pathsd
        keep = ones(Bool, length(paths))
        merged = false
        for i in 1:length(paths)
            for j in 1:(i-1)
                if keep[j] && canmerge(paths[i], paths[j], ztol)
                    #println("Merging path since ||x_i - x_j||_∞ = $(norm(getstatevalue(paths[j].sol) - getstatevalue(paths[i].sol), Inf))")
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

function childjobs(g::AbstractStochasticProgram, pathsd::Vector{Tuple{NodeT, Vector{SDDPPath{SolT}}}}, pathsampler, t, num_stages, endedpaths) where {SolT, NodeT}
    jobsd = Dict{NodeT, Vector{SDDPJob{SolT, NodeT}}}()
    for (node, paths) in pathsd
        if !isleaf(g, node)
            for path in paths
                # Adding Jobs
                npaths = samplepaths(pathsampler, g, node, path.K, t, num_stages)
                childocuts = Array{Any}(outdegree(g, node))
                for (i, child) in enumerate(out_neighbors(g, node))
                    if sum(npaths[i]) != 0 || needallchildsol(cutgenerator(g, node)) # || t == 2
                        addjob!(jobsd, child, SDDPJob{NodeT}(path.proba * probability(g, Edge(node, child)), npaths[i], node, path, child))
                    end
                end
            end
        else
            append!(endedpaths, paths)
        end
    end
    jobsd
end

function jobstopaths(jobsd::Dict{NodeT, Vector{SDDPJob{SolT, NodeT}}}, g::AbstractStochasticProgram) where {SolT, NodeT}
    pathsd = Tuple{NodeT, Vector{SDDPPath{SolT}}}[]
    for (node, jobs) in jobsd
        K = [find(job.K .!= 0) for job in jobs]
        keep = find(Bool[jobs[i].parent.childs_feasible && !isempty(K[i]) for i in 1:length(jobs)])
        if !isempty(keep)
            paths = SDDPPath{SolT}[SDDPPath(get(jobs[i].sol), jobs[i].parent.z[K[i]]+get(jobs[i].sol).objvalx, jobs[i].proba[K[i]], jobs[i].K[K[i]], outdegree(g, node)) for i in keep]
            push!(pathsd, (node, paths))
        end
    end
    pathsd
end

function solvejob!(sp::AbstractStochasticProgram, job::SDDPJob, node, stats)
    stats.setxtime += @_time setchildx!(sp, job.parentnode, job.child, job.parent.sol)
    stats.nsetx += 1
    stats.solvertime += @_time job.sol = solve!(sp, node)
    job.parent.childsols[job.child] = get(job.sol)
    stats.nsolved += 1
    if get(job.sol).status == :Infeasible
        job.parent.childs_feasible = false
    elseif get(job.sol).status == :Unbounded
        job.parent.childs_bounded = false
    end
end
