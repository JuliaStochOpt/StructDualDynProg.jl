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
    σ = sqrt(dot(proba, (z .- μ).^2))
    μ, σ
end

mutable struct SolutionPool{TT<:SOI.AbstractTransition, SolT<:SOI.AbstractSolution} <: SOI.AbstractSolutionPool
    # Parent solution
    parent_solution::SolT
    # Are all the child jobs solved feasible ?
    children_feasible::Bool
    # Are all the child jobs solved bounded ?
    children_bounded::Bool
    # List of children solutions
    children_solutions::Dict{TT, SolT}
end
SolutionPool{TT}(sol) where TT = SolutionPool{TT, typeof(sol)}(sol, true, true, Dict{TT, typeof(sol)}())
SOI.getsolution(pool::SolutionPool) = pool.parent_solution
SOI.allfeasible(pool::SolutionPool) = pool.children_feasible
SOI.allbounded(pool::SolutionPool) = pool.children_bounded
SOI.hassolution(pool::SolutionPool{TT}, tr::TT) where TT = haskey(pool.children_solutions, tr)
SOI.getsolution(pool::SolutionPool{TT}, tr::TT) where TT = pool.children_solutions[tr]

mutable struct SDDPPath{TT<:SOI.AbstractTransition, SolT<:SOI.AbstractSolution}
    # Solution pool for the last node of the path
    pool::SolutionPool{TT, SolT}

    # The three following arguments are vectors since when we merge paths,
    # we need to keep track of the difference of `z` that they have accumulated
    # while diverging for computing the variance `σ^2` in [`meanstdpaths`](@ref) at the end.
    # Sum of the objective solution (without θ) at each node of the path
    z::Vector{Float64}
    # Probability of the path
    proba::Vector{Float64}
    # Number of walks in the path
    K::Vector{Int}

    function SDDPPath{TT}(sol::SolT, z, proba, K) where {TT, SolT}
        new{TT, SolT}(SolutionPool{TT}(sol), z, proba, K)
    end
end

function meanstdpaths(paths::Vector{<:SDDPPath}, Ktot)
    z = Compat.reduce(append!, Vector{Float64}[x.z for x in paths], init=Float64[])
    proba = Compat.reduce(append!, Vector{Float64}[x.proba for x in paths], init=Float64[])
    npaths = Compat.reduce(append!, Vector{Int}[x.K for x in paths], init=Int[])
    meanstdpaths(z, proba, npaths, Ktot)
end

function Base.isapprox(p::SDDPPath, q::SDDPPath)
    Base.isapprox(SOI.getnodevalue(SOI.getsolution(p.pool)), SOI.getnodevalue(SOI.getsolution(q.pool)))
end

function canmerge(p::SDDPPath, q::SDDPPath, ztol)
    _isapprox(SOI.getnodevalue(SOI.getsolution(p.pool)), SOI.getnodevalue(SOI.getsolution(q.pool)), ztol)
end

function merge!(p::SDDPPath, q::SDDPPath)
    @assert p.pool.children_feasible == q.pool.children_feasible
    append!(p.z, q.z)
    append!(p.proba, q.proba)
    append!(p.K, q.K)
end

mutable struct Job{SolT<:SOI.AbstractSolution, NodeT, TT<:SOI.AbstractTransition}
    sol::Union{Nothing, SolT}
    proba::Vector{Float64}
    K::Vector{Int}
    parentnode::NodeT
    parent::SDDPPath{TT, SolT}
    tr::TT

    function Job(proba::Vector{Float64}, K::Vector{Int}, parentnode::NodeT, parent::SDDPPath{TT, SolT}, tr::TT) where {SolT, NodeT, TT}
        new{SolT, NodeT, TT}(nothing, proba, K, parentnode, parent, tr)
    end
end

function addjob!(jobsd::Dict{NodeT, Vector{Job{SolT, NodeT, TT}}}, node::NodeT, job::Job{SolT, NodeT, TT}) where {SolT, NodeT, TT}
    if node in keys(jobsd)
        push!(jobsd[node], job)
    else
        jobsd[node] = [job]
    end
end

"""
    mergesamepaths(pathsd::Vector{Tuple{NodeT, Vector{<:SDDPPath}}}, to::TimerOutput, ztol) where NodeT

Find paths that are at the same node with the same parent solution and merge them.
It could happend that two path diverges (go to different node) but then meet again at the same node, if their parent solution is the same then each path will do exactly the same computation so merging them should save time.
"""
function mergesamepaths(pathsd::Vector{Tuple{NodeT, Vector{SDDPPath{TT, SolT}}}}, to::TimerOutput, ztol) where {SolT, NodeT, TT}
    before = sum([sum([sum(path.K) for path in paths]) for (node, paths) in pathsd])
    newpathsd = Tuple{NodeT, Vector{SDDPPath{TT, SolT}}}[]
    @timeit to "merged" for (node, paths) in pathsd
        keep = ones(Bool, length(paths))
        merged = false
        for i in 1:length(paths)
            for j in 1:(i-1)
                if keep[j] && canmerge(paths[i], paths[j], ztol)
                    #println("Merging path since ||x_i - x_j||_∞ = $(norm(getnodevalue(paths[j].sol) - getnodevalue(paths[i].sol), Inf))")
                    merge!(paths[i], paths[j])
                    keep[j] = false
                    merged = true
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
    childjobs(g::SOI.AbstractStochasticProgram, pathsd::Vector{Tuple{NodeT, Vector{<:SDDPPath}}}, pathsampler::AbstractPathSampler, t, num_stages, endedpaths) where SolT

Given paths in `pathsd`, put the paths that have no child in `endedpaths` and sample child jobs using `pathsample` for other paths.
"""
function childjobs(sp::SOI.AbstractStochasticProgram, pathsd::Vector{Tuple{NodeT, Vector{SDDPPath{TT, SolT}}}}, pathsampler::AbstractPathSampler, t, num_stages, endedpaths) where {SolT, NodeT, TT}
    jobsd = Dict{NodeT, Vector{Job{SolT, NodeT, SOI.get(sp, SOI.TransitionType())}}}()
    for (node, paths) in pathsd
        if !isempty(SOI.get(sp, SOI.OutTransitions(), node))
            for path in paths
                # Adding Jobs
                npaths = samplepaths(pathsampler, sp, node, path.K, t, num_stages)
                for (i, tr) in enumerate(SOI.get(sp, SOI.OutTransitions(), node))
                    if !iszero(sum(npaths[i])) || SOI.get(sp, SOI.NeedAllSolutions(), node) # || t == 2
                        addjob!(jobsd, SOI.get(sp, SOI.Target(), tr), Job(path.proba * SOI.get(sp, SOI.Probability(), tr), npaths[i], node, path, tr))
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
    jobstopath(jobsd::Dict{NodeT, Vector{<:Job}}, g::SOI.AbstractStochasticProgram)

Transforms the jobs `jobsd` created by [`childjobs`](@ref) to to paths.
"""
function jobstopaths(jobsd::Dict{NodeT, Vector{Job{SolT, NodeT, TT}}}, g::SOI.AbstractStochasticProgram) where {SolT, NodeT, TT}
    pathsd = Tuple{NodeT, Vector{SDDPPath{TT, SolT}}}[]
    for (node, jobs) in jobsd
        # We create a job even if there is no path going to the node in case
        # we want to create an AveragedCut (since in this case we need to solve all children).
        # However we do not want to create a path for these jobs so we filter them out.
        K = [findall(job.K .!= 0) for job in jobs]
        keep = findall(Bool[jobs[i].parent.pool.children_feasible && !isempty(K[i]) for i in 1:length(jobs)])
        if !isempty(keep)
            paths = SDDPPath{TT, SolT}[SDDPPath{TT}(jobs[i].sol, jobs[i].parent.z[K[i]] .+ jobs[i].sol.objvalx, jobs[i].proba[K[i]], jobs[i].K[K[i]]) for i in keep]
            push!(pathsd, (node, paths))
        end
    end
    pathsd
end

"""
    solvejob!(sp::SOI.AbstractStochasticProgram, job::Job, node, to::TimerOutput)

Solves the job `job` of node `node`.
"""
function solvejob!(sp::SOI.AbstractStochasticProgram, job::Job, node, to::TimerOutput)
    @timeit to "setx" SOI.set!(sp, SOI.SourceSolution(), job.tr, SOI.getsolution(job.parent.pool))
    @timeit to "solve" job.sol = SOI.get(sp, SOI.Solution(), node)
    job.parent.pool.children_solutions[job.tr] = job.sol
    if job.sol.status == :Infeasible
        job.parent.pool.children_feasible = false
    elseif job.sol.status == :Unbounded
        job.parent.pool.children_bounded = false
    end
end
