type SDDPSolution
    status
    objval
    sol
    attrs
end

macro mytime(x)
    quote
        y = @timed $(esc(x))
        # y[1] is returned value
        # y[2] is time in seconds
        y[2]
    end
end

type SDDPStats
    nsolved::Int
    solvertime::Float64
    nmerged::Int
    mergetime::Float64
    nfcuts::Int
    fcutstime::Float64
    nocuts::Int
    ocutstime::Float64
    nsetx::Int
    setxtime::Float64
end

SDDPStats() = SDDPStats(0,.0,0,.0,0,.0,0,.0,0,.0)

import Base: +, show

function +(a::SDDPStats, b::SDDPStats)
    SDDPStats(a.nsolved + b.nsolved, a.solvertime + b.solvertime,
    a.nmerged + b.nmerged, a.mergetime  + b.mergetime,
    a.nfcuts  + b.nfcuts , a.fcutstime  + b.fcutstime,
    a.nocuts  + b.nocuts , a.ocutstime  + b.ocutstime,
    a.nsetx   + b.nsetx  , a.setxtime   + b.setxtime)
end

function showtime(t::Float64)
    if !isfinite(t)
        "   min    s    ms    μs"
    else
        s = Int(floor(t))
        t = (t - s) * 1000
        min = div(s, 60)
        s = mod(s, 60)
        ms = Int(floor(t))
        t = (t - ms) * 1000
        μs = Int(floor(t))
        @sprintf "%3dmin %3ds %3dms %3dμs" min s ms μs
    end
end

function Base.show(io::IO, stat::SDDPStats)
    println("                        |     Total time [s]      |  Number  | Average time [s]")
    @printf "        Solving problem | %s | %8d | %s\n" showtime(stat.solvertime) stat.nsolved showtime(stat.solvertime / stat.nsolved)
    @printf "          Merging paths | %s | %8d | %s\n" showtime(stat.mergetime ) stat.nmerged showtime(stat.mergetime  / stat.nmerged)
    @printf "Adding feasibility cuts | %s | %8d | %s\n" showtime(stat.fcutstime ) stat.nfcuts  showtime(stat.fcutstime  / stat.nfcuts )
    @printf "Adding  optimality cuts | %s | %8d | %s\n" showtime(stat.ocutstime ) stat.nocuts  showtime(stat.ocutstime  / stat.nocuts )
    @printf "Setting parent solution | %s | %8d | %s\n" showtime(stat.setxtime  ) stat.nsetx   showtime(stat.setxtime   / stat.nsetx  )
    @printf "                        | %s |" showtime(stat.solvertime+stat.fcutstime+stat.ocutstime+stat.setxtime)
end

function meanstdpaths(z::Vector{Float64}, proba::Vector{Float64}, npaths::Vector{Int}, Ktot)
    if Ktot != -1
        @assert sum(npaths) == Ktot
        proba = npaths / Ktot
    end
    μ = dot(proba, z)
    σ = sqrt(dot(proba, (z - μ).^2))
    μ, σ
end

function choosepaths(node::SDDPNode, K::Int, pathsel, t, num_stages)
    if K == -1
        map(child->-1, node.children)
    else
        if pathsel == :nPaths
            den = numberofpaths(node, t-1, num_stages)
            pmf = map(child->numberofpaths(child, t, num_stages) / den, node.children)
        else
            pmf = node.proba
        end
        cmf = cumsum(pmf)
        @assert abs(cmf[end] - 1) < 1e-6
        cmf[end] = 1
        samples = rand(Float64, K)
        npaths = zeros(Int, length(node.children))
        sort!(samples)
        i = 1
        for j in samples
            while j >= cmf[i]
                i += 1
            end
            npaths[i] += 1
        end
        npaths
    end
end

function choosepaths(node::SDDPNode, K::Vector{Int}, pathsel, t, num_stages)
    npathss = Vector{Int}[similar(K) for i in 1:length(node.children)]
    for i in 1:length(K)
        npaths = choosepaths(node, K[i], pathsel, t, num_stages)
        for c in 1:length(node.children)
            npathss[c][i] = npaths[c]
        end
    end
    npathss
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
    myeq(p.sol.x, q.sol.x, ztol)
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

function addjob!(jobsd::Dict{SDDPNode, Vector{SDDPJob}}, node::SDDPNode, job::SDDPJob)
    if node in keys(jobsd)
        push!(jobsd[node], job)
    else
        jobsd[node] = [job]
    end
end

function iteration{S}(root::SDDPNode{S}, Ktot::Int, num_stages, verbose, pathsel, ztol)
    stats = SDDPStats()

    stats.solvertime += @mytime rootsol = loadAndSolve(root)
    stats.nsolved += 1
    infeasibility_detected = rootsol.status == :Infeasible
    if infeasibility_detected
        pathsd = Dict{SDDPNode, Vector{SDDPPath}}()
    else
        pathsd = Dict{SDDPNode, Vector{SDDPPath}}(root => [SDDPPath(rootsol, [rootsol.objvalx], [1.], [Ktot], length(root.children))])
    end
    endedpaths = SDDPPath[]

    for t in 2:num_stages
        if verbose >= 3
            @show t
        end

        # Merge paths
        if true
            before = sum([sum([sum(path.K) for path in paths]) for (node, paths) in pathsd])
            stats.mergetime += @mytime for (node, paths) in pathsd
                keep = ones(Bool, length(paths))
                merged = false
                for i in 1:length(paths)
                    for j in 1:(i-1)
                        if keep[j] && canmerge(paths[i], paths[j], ztol)
                            @show maximum(abs(paths[j].sol.x - paths[i].sol.x))
                            merge!(paths[i], paths[j])
                            keep[j] = false
                            merged = true
                            stats.nmerged += 1
                            break
                        end
                    end
                end
                pathsd[node] = paths[keep]
            end
            after = sum([sum([sum(path.K) for path in paths]) for (node, paths) in pathsd])
            @assert before == after
        end

        # Make jobs
        jobsd = Dict{SDDPNode, Vector{SDDPJob}}()
        for (parent, paths) in pathsd
            if isempty(parent.children)
                append!(endedpaths, paths)
            else
                for path in paths
                    # Adding Jobs
                    npaths = choosepaths(parent, path.K, pathsel, t, num_stages)
                    childocuts = Array{Any}(length(parent.children))
                    for i in 1:length(parent.children)
                        if t == 2 || sum(npaths[i]) != 0 || parent.nlds.cutmode == :AveragedCut
                            addjob!(jobsd, parent.children[i], SDDPJob(path.proba * parent.proba[i], npaths[i], parent, path, i))
                        end
                    end
                end
            end
        end

        # Solve Jobs (parallelism possible here)
        for (node, jobs) in jobsd
            for job in jobs
                if job.parent.childs_feasible
                    stats.setxtime += @mytime setchildx(job.parentnode, job.i, job.parent.sol)
                    stats.nsetx += 1
                    stats.solvertime += @mytime job.sol = loadAndSolve(node)
                    job.parent.childsols[job.i] = job.sol
                    stats.nsolved += 1
                    if get(job.sol).status == :Infeasible
                        job.parent.childs_feasible = false
                    elseif get(job.sol).status == :Unbounded
                        job.parent.childs_bounded = false
                    end
                end
            end
        end

        # Add cuts
        # Feasibility cut
        # D = π T
        # d = π h + σ d
        # Optimality cut
        # E = π T
        # e = π h + ρ e + σ d
        for (parent, paths) in pathsd
            for path in paths
                if parent.nlds.cutmode == :AveragedCut && path.childs_feasible && path.childs_bounded
                    avga = zeros(parent.nlds.nx)
                    avgβ = 0
                end
                for i in 1:length(parent.children)
                    if isnull(path.childsols[i])
                        @assert !path.childs_feasible || parent.nlds.cutmode != :AveragedCut
                    elseif get(path.childsols[i]).status != :Unbounded
                        childsol = get(path.childsols[i])
                        a = childsol.πT
                        β = childsol.πh + childsol.σd
                        if !path.childs_feasible
                            if childsol.status == :Infeasible
                                infeasibility_detected = true
                                stats.nfcuts += 1
                                stats.fcutstime += @mytime pushfeasibilitycut!(parent.children[i], a, β, parent)
                            end
                        elseif !(parent.nlds.cutmode == :AveragedCut && (!path.childs_bounded || !path.childs_feasible))
                            @assert childsol.status == :Optimal
                            if isnull(parent.childT)
                                aT = a
                            else
                                aT = get(parent.childT)[i]' * a
                            end
                            β += childsol.ρe
                            # This assumption is dropped now
                           #if mylt(β - dot(aT, path.sol.x), .0, ztol)
                           #    error("The objectives are supposed to be nonnegative")
                           #end
                            if parent.nlds.cutmode == :MultiCut
                                if path.sol.status == :Unbounded || mylt(path.sol.θ[i], β - dot(aT, path.sol.x), ztol)
                                    stats.ocutstime += @mytime pushoptimalitycutforparent!(parent.children[i], a, β, parent)
                                    stats.nocuts += 1
                                end
                            elseif parent.nlds.cutmode == :AveragedCut
                                avga += parent.proba[i] * aT
                                avgβ += parent.proba[i] * β
                            end
                        end
                    end
                end
                if parent.nlds.cutmode == :AveragedCut && path.childs_feasible && path.childs_bounded
                    if path.sol.status == :Unbounded || mylt(path.sol.θ[1], avgβ - dot(avga, path.sol.x), ztol)
                        stats.ocutstime += @mytime pushoptimalitycut!(parent, avga, avgβ, parent)
                        stats.nocuts += 1
                    end
                end
            end
        end

        # Apply cut addition
        for (parent, paths) in pathsd
            for child in parent.children
                applyfeasibilitycut!(child)
            end
            if parent.nlds.cutmode == :MultiCut
                for child in parent.children
                    applyoptimalitycutforparent!(child)
                end
            elseif parent.nlds.cutmode == :AveragedCut
                applyoptimalitycut!(parent)
            end
        end

        # Jobs -> Paths
        empty!(pathsd)
        for (node, jobs) in jobsd
            K = [find(job.K .!= 0) for job in jobs]
            keep = find(Bool[jobs[i].parent.childs_feasible && !isempty(K[i]) for i in 1:length(jobs)])
            if !isempty(keep)
                paths = SDDPPath[SDDPPath(get(jobs[i].sol), jobs[i].parent.z[K[i]]+get(jobs[i].sol).objvalx, jobs[i].proba[K[i]], jobs[i].K[K[i]], length(node.children)) for i in keep]
                pathsd[node] = paths
            end
        end
    end

    if infeasibility_detected
        z_UB = Inf # FIXME assumes minimization
        σ = 0
    else
        for (node, paths) in pathsd
            append!(endedpaths, paths)
        end
        z_UB, σ = meanstdpaths(endedpaths, Ktot)
    end
    rootsol, stats, z_UB, σ
end

"""
$(SIGNATURES)

Runs the SDDP algorithms on the lattice given by `root`.
The algorithm will do iterations until `stopcrit` decides to stop or when the root node is infeasible.
In each iterations, `K` paths will be explored up to `num_stages` stages.
The paths will be selected according to `pathsel` and equivalent paths might be merged if their difference is smaller than `ztol`.
The parameter `ztol` is also used to check whether a new cut is useful.
"""
function SDDP(root::SDDPNode, num_stages; K::Int=25, stopcrit::AbstractStoppingCriterion=Pereira(), verbose=0, pathsel::Symbol=:Proba, ztol=1e-6)
    if !(pathsel in [:Proba, :nPaths])
        error("Invalid pathsel")
    end
    rootsol = nothing
    totaltime = 0
    totalstats = SDDPStats()

    z_UB = Inf
    σ = 0
    z_LB = 0
    iter = 0
    nfcuts = 0
    nocuts = 0

    while (rootsol === nothing || rootsol.status != :Infeasible) && !stop(stopcrit, iter, nfcuts, nocuts, K, z_LB, z_UB, σ)
        itertime = @mytime rootsol, stats, z_UB, σ = iteration(root, K, num_stages, verbose, pathsel, ztol)
        z_LB = rootsol.objval
        iter += 1
        nfcuts = stats.nfcuts
        nocuts = stats.nocuts

        totaltime  += itertime
        totalstats += stats
        if verbose >= 2
            println("Iteration $iter completed in $itertime s (Total time is $totaltime)")
            println("Status: $(rootsol.status)")
            println("z_UB: $(z_UB)")
            println("z_LB: $(z_LB)")
            #println(" Solution value: $(rootsol.x)")
            println("Stats for this iteration:")
            println(stats)
            println("Total stats:")
            println(totalstats)
        end
    end

    if verbose >= 1
        println("SDDP completed in $iter iterations in $totaltime s")
        println("Status: $(rootsol.status)")
        #println("Objective value: $(rootsol.objval)")
        println("z_UB: $(z_UB)")
        println("z_LB: $(z_LB)")
        #println(" Solution value: $(rootsol.x)")
        println("Total stats:")
        println(totalstats)
    end

    attrs = Dict()
    attrs[:niter] = iter
    attrs[:nfcuts] = nfcuts
    attrs[:nocuts] = nocuts
    SDDPSolution(rootsol.status, rootsol.objval, rootsol.x, attrs)
end
