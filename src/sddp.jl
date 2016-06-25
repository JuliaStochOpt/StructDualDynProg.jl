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

function meanstdpaths(z::Vector{Float64}, proba::Vector{Float64}, npaths::Vector{Int}, totalmccount)
  if totalmccount != -1
    @assert sum(npaths) == totalmccount
    proba = npaths / totalmccount
  end
  μ = dot(proba, z)
  σ = sqrt(dot(proba, (z - μ).^2))
  μ, σ
end

function choosepaths(node::SDDPNode, mccount::Int, pathsel, t, num_stages)
  if mccount == -1
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
    samples = rand(Float64, mccount)
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

function choosepaths(node::SDDPNode, mccount::Vector{Int}, pathsel, t, num_stages)
  npathss = Vector{Int}[similar(mccount) for i in 1:length(node.children)]
  for i in 1:length(mccount)
    npaths = choosepaths(node, mccount[i], pathsel, t, num_stages)
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
  mccount::Vector{Int}
  feasible::Bool
  childsols::Vector{Nullable{NLDSSolution}}

  function SDDPPath(sol, z, proba, mccount, nchilds)
    childsols = Nullable{NLDSSolution}[nothing for i in 1:nchilds]
    new(sol, z, proba, mccount, true, childsols)
  end
end

function meanstdpaths(paths::Vector{SDDPPath}, totalmccount)
  z = reduce(append!, Float64[], Vector{Float64}[x.z for x in paths])
  proba = reduce(append!, Float64[], Vector{Float64}[x.proba for x in paths])
  npaths = reduce(append!, Int[], Vector{Int}[x.mccount for x in paths])
  meanstdpaths(z, proba, npaths, totalmccount)
end

function canmerge(p::SDDPPath, q::SDDPPath, ztol)
  myeq(p.sol.x, q.sol.x, ztol)
end

function merge!(p::SDDPPath, q::SDDPPath)
  @assert p.feasible == q.feasible
  append!(p.z, q.z)
  append!(p.proba, q.proba)
  append!(p.mccount, q.mccount)
end

type SDDPJob
  sol::Nullable{NLDSSolution}
  proba::Vector{Float64}
  mccount::Vector{Int}
  parentnode::SDDPNode
  parent::SDDPPath
  i::Int

  function SDDPJob(proba::Vector{Float64}, mccount::Vector{Int}, parentnode::SDDPNode, parent::SDDPPath, i::Int)
    new(nothing, proba, mccount, parentnode, parent, i::Int)
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

function iteration(root::SDDPNode, totalmccount::Int, num_stages, verbose, pathsel, ztol)
  stats = SDDPStats()

  stats.solvertime += @mytime rootsol = loadAndSolve(root)
  stats.nsolved += 1
  infeasibility_detected = rootsol.status == :Infeasible
  if infeasibility_detected
    pathsd = Dict{SDDPNode, Vector{SDDPPath}}()
  else
    pathsd = Dict{SDDPNode, Vector{SDDPPath}}(root => [SDDPPath(rootsol, [rootsol.objvalx], [1.], [totalmccount], length(root.children))])
  end
  endedpaths = SDDPPath[]

  for t in 2:num_stages
    if verbose >= 3
      @show t
    end

    # Merge paths
    if true
      before = sum([sum([sum(path.mccount) for path in paths]) for (node, paths) in pathsd])
      stats.mergetime += @mytime for (node, paths) in pathsd
        keep = ones(Bool, length(paths))
        merged = false
        for i in 1:length(paths)
          for j in 1:(i-1)
            if keep[j] && canmerge(paths[i], paths[j], ztol)
              #@show paths[i].sol.x
              #@show paths[j].sol.x
              #@show paths[j].sol.x - paths[i].sol.x
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
      after = sum([sum([sum(path.mccount) for path in paths]) for (node, paths) in pathsd])
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
          npaths = choosepaths(parent, path.mccount, pathsel, t, num_stages)
          childocuts = Array{Any}(length(parent.children))
          for i in 1:length(parent.children)
            if t == 2 || sum(npaths[i]) != 0 || parent.nlds.cutmode == :AveragedCut
              addjob!(jobsd, parent.children[i], SDDPJob(path.proba * parent.proba[i], npaths[i], parent, path, i))
            end
          end
        end
      end
    end

    # Solving Jobs (paralellism possible here)
    for (node, jobs) in jobsd
      for job in jobs
        if job.parent.feasible
          stats.setxtime += @mytime setchildx(job.parentnode, job.i, job.parent.sol.x)
          stats.nsetx += 1
          stats.solvertime += @mytime job.sol = loadAndSolve(node)
          job.parent.childsols[job.i] = job.sol
          stats.nsolved += 1
          if get(job.sol).status == :Infeasible
            job.parent.feasible = false
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
        if parent.nlds.cutmode == :AveragedCut && path.feasible
          avga = zeros(parent.nlds.nx)
          avgβ = 0
        end
        for i in 1:length(parent.children)
          if isnull(path.childsols[i])
            @assert !path.feasible || parent.nlds.cutmode != :AveragedCut
          else
            childsol = get(path.childsols[i])
            a = childsol.πT
            β = childsol.πh + childsol.σd
            if !path.feasible
              if childsol.status == :Infeasible
                infeasibility_detected = true
                stats.nfcuts += 1
                stats.fcutstime += @mytime pushfeasibilitycut!(parent.children[i], a, β, parent)
              end
            else
              @assert childsol.status == :Optimal
              if isnull(parent.childT)
                aT = a
              else
                aT = get(parent.childT)[i]' * a
              end
              β += childsol.ρe
              if mylt(β - dot(aT, path.sol.x), .0, ztol)
                error("The objectives are supposed to be nonnegative")
              end
              if parent.nlds.cutmode == :MultiCut
                if mylt(path.sol.θ[i], β - dot(aT, path.sol.x), ztol)
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
        if parent.nlds.cutmode == :AveragedCut && path.feasible
          if mylt(path.sol.θ[1], avgβ - dot(avga, path.sol.x), ztol)
            stats.ocutstime += @mytime pushoptimalitycut!(parent, avga, avgβ, parent)
            stats.nocuts += 1
          end
        end
      end
    end

    # Jobs -> Paths
    empty!(pathsd)
    for (node, jobs) in jobsd
      K = [find(job.mccount .!= 0) for job in jobs]
      keep = find(Bool[jobs[i].parent.feasible && !isempty(K[i]) for i in 1:length(jobs)])
      if !isempty(keep)
        paths = SDDPPath[SDDPPath(get(jobs[i].sol), jobs[i].parent.z[K[i]]+get(jobs[i].sol).objvalx, jobs[i].proba[K[i]], jobs[i].mccount[K[i]], length(node.children)) for i in keep]
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
    z_UB, σ = meanstdpaths(endedpaths, totalmccount)
  end
  rootsol, stats, z_UB, σ
end

function SDDP(root::SDDPNode, num_stages, mccount::Int=25, verbose=0, pereiracoef=2, stopcrit::Function=(x,y)->false, pathsel::Symbol=:Proba, ztol=1e-10)
  if !(pathsel in [:Proba, :nPaths])
    error("Invalid pathsel")
  end
  cut_added = true
  niter = 0
  nfcuts = 0
  nocuts = 0
  rootsol = nothing
  totaltime = 0
  totalstats = SDDPStats()

  z_LB = 0
  z_UB = Inf
  σ = 0
  stdmccoef = 0.05

  while (mccount != -1 || cut_added) && (mccount == -1 || z_LB < z_UB - pereiracoef * σ / sqrt(mccount) || σ / sqrt(mccount) > stdmccoef * z_LB) && (rootsol === nothing || rootsol.status != :Infeasible) && (niter == 0 || !stopcrit(niter, rootsol.objval))
    niter += 1
    cut_added = false
    itertime = @mytime rootsol, stats, z_UB, σ = iteration(root, mccount, num_stages, verbose, pathsel, ztol)
    #Lower bound since θ >= 0
    z_LB = rootsol.objval

    totaltime  += itertime
    totalstats += stats
    if verbose >= 2
      println("Iteration $niter completed in $itertime s (Total time is $totaltime)")
      println("Status: $(rootsol.status)")
      println("z_UB: $(z_UB)")
      println("z_LB: $(z_LB)")
      if mccount != -1
        println("pereira bound: $(z_UB - pereiracoef * σ / sqrt(mccount))")
      end
      #println(" Solution value: $(rootsol.x)")
      println("Stats for this iteration:")
      println(stats)
      println("Total stats:")
      println(totalstats)
    end
    cut_added = stats.nfcuts > 0 || stats.nocuts > 0
  end

  if verbose >= 1
    println("SDDP completed in $niter iterations in $totaltime s")
    println("Status: $(rootsol.status)")
    #println("Objective value: $(rootsol.objval)")
    println("z_UB: $(z_UB)")
    println("z_LB: $(z_LB)")
    if mccount != -1
      println("pereira bound: $(z_UB - pereiracoef * σ / sqrt(mccount))")
    end
    #println(" Solution value: $(rootsol.x)")
    println("Total stats:")
    println(totalstats)
  end

  attrs = Dict()
  attrs[:niter] = niter
  attrs[:nfcuts] = nfcuts
  attrs[:nocuts] = nocuts
  SDDPSolution(rootsol.status, rootsol.objval, rootsol.x, attrs)
end
