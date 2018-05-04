function getproba(proba::Vector{Float64}, npaths::Vector{Int}, Ktot)
    if Ktot == -1
        # We have thrown one walk per possible path so npaths does not represent the probability.
        # Therefore we use proba instead.
        proba
    else
        # npaths is representative since the paths were sampled using the probabilities.
        @assert sum(npaths) == Ktot
        npaths ./ Ktot
    end
end
function meanstdpaths(z::Vector{Float64}, proba::Vector{Float64}, npaths::Vector{Int}, Ktot)
    proba = getproba(proba, npaths, Ktot)
    μ = dot(proba, z)
    σ = sqrt(dot(proba, (z - μ).^2))
    μ, σ
end

mutable struct SDDPPath{SolT}
    # Solution of the last node of the path
    sol::SolT
    # The three following arguments are vectors since when we merge paths,
    # we need to keep track of the difference of `z` that they have accumulated
    # while diverging for computing the variance `σ^2` in [`meanstdpaths`](@ref) at the end.
    # Sum of the objective solution (without θ) at each node of the path
    z::Vector{Float64}
    # Probability of the path
    proba::Vector{Float64}
    # Number of walks in the path
    K::Vector{Int}
    # Are all the child jobs solved feasible ?
    childs_feasible::Bool
    # Are all the child jobs solved bounded ?
    childs_bounded::Bool
    # List of children solutions
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

mutable struct Job{SolT, NodeT}
    sol::Nullable{SolT}
    proba::Vector{Float64}
    K::Vector{Int}
    parentnode::NodeT
    parent::SDDPPath{SolT}
    child::NodeT

    function Job{NodeT}(proba::Vector{Float64}, K::Vector{Int}, parentnode::NodeT, parent::SDDPPath{SolT}, child::NodeT) where {SolT, NodeT}
        new{SolT, NodeT}(nothing, proba, K, parentnode, parent, child)
    end
end

function addjob!(jobsd::Dict{NodeT, Vector{Job{SolT, NodeT}}}, node::NodeT, job::Job{SolT, NodeT}) where {SolT, NodeT}
    if node in keys(jobsd)
        push!(jobsd[node], job)
    else
        jobsd[node] = [job]
    end
end

"""
    mergesamepaths(pathsd::Vector{Tuple{NodeT, Vector{SDDPPath{SolT}}}}, stats, ztol) where {SolT, NodeT}

Find paths that are at the same node with the same parent solution and merge them.
It could happend that two path diverges (go to different node) but then meet again at the same node, if their parent solution is the same then each path will do exactly the same computation so merging them should save time.
"""
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

"""
    childjobs(g::AbstractStochasticProgram, pathsd::Vector{Tuple{NodeT, Vector{SDDPPath{SolT}}}}, pathsampler, t, num_stages, endedpaths) where {SolT, NodeT}

Given paths in `pathsd`, put the paths that have no child in `endedpaths` and sample child jobs using `pathsample` for other paths.
"""
function childjobs(g::AbstractStochasticProgram, pathsd::Vector{Tuple{NodeT, Vector{SDDPPath{SolT}}}}, pathsampler, t, num_stages, endedpaths) where {SolT, NodeT}
    jobsd = Dict{NodeT, Vector{Job{SolT, NodeT}}}()
    for (node, paths) in pathsd
        if !isleaf(g, node)
            for path in paths
                # Adding Jobs
                npaths = samplepaths(pathsampler, g, node, path.K, t, num_stages)
                for (i, tr) in enumerate(out_transitions(g, node))
                    child = target(tr)
                    if !iszero(sum(npaths[i])) || needallchildsol(cutgenerator(g, node)) # || t == 2
                        addjob!(jobsd, child, Job{NodeT}(path.proba * probability(g, Edge(node, child)), npaths[i], node, path, child))
                    end
                end
            end
        else
            append!(endedpaths, paths)
        end
    end
    jobsd
end

"""
    jobstopath(jobsd::Dict{NodeT, Vector{Job{SolT, NodeT}}}, g::AbstractStochasticProgram) where {SolT, NodeT}

Transforms the jobs `jobsd` created by [`childjobs`](@ref) to to paths.
"""
function jobstopaths(jobsd::Dict{NodeT, Vector{Job{SolT, NodeT}}}, g::AbstractStochasticProgram) where {SolT, NodeT}
    pathsd = Tuple{NodeT, Vector{SDDPPath{SolT}}}[]
    for (node, jobs) in jobsd
        # We create a job even if there is no path going to the node in case
        # we want to create an AveragedCut (since in this case we need to solve all children).
        # However we do not want to create a path for these jobs so we filter them out.
        K = [find(job.K .!= 0) for job in jobs]
        keep = find(Bool[jobs[i].parent.childs_feasible && !isempty(K[i]) for i in 1:length(jobs)])
        if !isempty(keep)
            paths = SDDPPath{SolT}[SDDPPath(get(jobs[i].sol), jobs[i].parent.z[K[i]]+get(jobs[i].sol).objvalx, jobs[i].proba[K[i]], jobs[i].K[K[i]], outdegree(g, node)) for i in keep]
            push!(pathsd, (node, paths))
        end
    end
    pathsd
end

"""
    solvejob!(sp::AbstractStochasticProgram, job::Job, node, stats)

Solves the job `job` of node `node`.
"""
function solvejob!(sp::AbstractStochasticProgram, job::Job, node, stats)
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
