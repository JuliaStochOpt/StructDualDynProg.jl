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
  nfcuts::Int
  fcutstime::Float64
  nocuts::Int
  ocutstime::Float64
  nsetx::Int
  setxtime::Float64
end

SDDPStats() = SDDPStats(0,.0,0,.0,0,.0,0,.0)

import Base: +, show

function +(a::SDDPStats, b::SDDPStats)
  SDDPStats(a.nsolved + b.nsolved, a.solvertime + b.solvertime,
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
  println("                        |     Total time [s]      | Number | Average time [s]")
  @printf "        Solving problem | %s | %6d | %s\n" showtime(stat.solvertime) stat.nsolved showtime(stat.solvertime / stat.nsolved)
  @printf "Adding feasibility cuts | %s | %6d | %s\n" showtime(stat.fcutstime ) stat.nfcuts  showtime(stat.fcutstime  / stat.nfcuts )
  @printf "Adding  optimality cuts | %s | %6d | %s\n" showtime(stat.ocutstime ) stat.nocuts  showtime(stat.ocutstime  / stat.nocuts )
  @printf "Setting parent solution | %s | %6d | %s\n" showtime(stat.setxtime  ) stat.nsetx   showtime(stat.setxtime   / stat.nsetx  )
  @printf "                        | %s |" showtime(stat.solvertime+stat.fcutstime+stat.ocutstime+stat.setxtime)
end

function choosepaths(node::SDDPNode, mccount, pathsel, t, num_stages)
  if mccount == :All
    map(child->:All, node.children)
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

function iteration{S}(root::SDDPNode{S}, totalmccount, num_stages, verbose, pathsel, TOL)
  rootsol = loadAndSolve(root)
  pathss = [(root, rootsol, rootsol.objvalx, 1., totalmccount)]

  stats = SDDPStats()
  infeasibility_detected = false
  for t in 2:num_stages
    if verbose >= 3
      @show t
    end
    # children are at t, parents are at t-1
    newpathss = []
    for (parent, psol, z, prob, mccount) in pathss
      npaths = choosepaths(parent, mccount, pathsel, t, num_stages)
      curpathss = [] # FIXME could already be allocated looking the number of nonzeros elements of npaths
      childsolved = zeros(Bool, length(parent.children))
      feasible = true
      nnewfcuts = 0
      nnewocuts = 0
      childocuts = Array{Any}(length(parent.children))
      for i in 1:length(parent.children)
        child = parent.children[i]
        if t == 2 || npaths[i] == :All || npaths[i] > 0 || parent.nlds.cutmode == :AveragedCut
          stats.setxtime += @mytime setchildx(parent, i, psol.x)
          stats.nsetx += 1
          stats.solvertime += @mytime childsol = loadAndSolve(child)
          stats.nsolved += 1
          childsolved[i] = true

          # Feasibility cut
          # D = π T
          # d = π h + σ d
          # Optimality cut
          # E = π T
          # e = π h + ρ e + σ d
          coef = childsol.πT
          rhs = childsol.πh + childsol.σd
          if childsol.status == :Infeasible
            infeasibility_detected = true
            feasible = false
            nnewfcuts += 1
            stats.fcutstime += @mytime pushfeasibilitycut!(child, coef, rhs, parent)
            break
          else
            rhs += childsol.ρe
            childocuts[i] = (coef, rhs)
          end
          if npaths[i] == :All || npaths[i] > 0
            push!(curpathss, (child, childsol, z + childsol.objvalx, prob * parent.proba[i], npaths[i]))
          end
        end
      end
      if feasible
        if parent.nlds.cutmode == :MultiCut
          for i in 1:length(parent.children)
            if childsolved[i]
              a, β = childocuts[i][1], childocuts[i][2]
              if psol.θ[i] < (β - dot(a, psol.x)) - TOL
                stats.ocutstime += @mytime pushoptimalitycutforparent!(parent.children[i], a, β, parent)
                nnewocuts += 1
              end
            end
          end
        elseif parent.nlds.cutmode == :AveragedCut
          if !isempty(parent.children)
            if isnull(parent.childT)
              a = sum(map(i->childocuts[i][1]*parent.proba[i], 1:length(parent.children)))
            else
              a = sum(map(i->get(parent.childT)[i]'*childocuts[i][1]*parent.proba[i], 1:length(parent.children)))
            end
            β = sum(map(i->childocuts[i][2]*parent.proba[i], 1:length(parent.children)))
            if psol.θ[1] < (β - dot(a, psol.x)) - TOL
              stats.ocutstime += @mytime pushoptimalitycut!(parent, a, β, parent)
              nnewocuts += 1
            end
          end
        end
      end
      stats.nfcuts += nnewfcuts
      stats.nocuts += nnewocuts
      if feasible && !isempty(curpathss)
        @assert nnewfcuts == 0
        append!(newpathss, curpathss)
      end
    end
    pathss = newpathss
  end
  if infeasibility_detected
    z_UB = Inf # FIXME assumes minimization
    σ = 0
  else
    z = Vector{Float64}(map(x->x[3], pathss))
    if totalmccount == :All
      prob = Vector{Float64}(map(x->x[4], pathss))
    else
      npathss = Vector{Int}(map(x->x[5], pathss))
      @assert sum(npathss) == totalmccount
      prob = npathss / totalmccount
    end
    z_UB = dot(prob, z)
    σ = sqrt(dot(prob, (z - z_UB).^2))
  end
  rootsol, stats, z_UB, σ
end

function SDDP(root::SDDPNode, num_stages, mccount=25, verbose=0, pereiracoef=2, stopcrit::Function=(x,y)->false, pathsel::Symbol=:Proba, TOL=1e-5)
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

  while (mccount != :All || cut_added) && (mccount == :All || z_LB < z_UB - pereiracoef * σ / sqrt(mccount) || σ / sqrt(mccount) > stdmccoef * z_LB) && (rootsol === nothing || rootsol.status != :Infeasible) && (niter == 0 || !stopcrit(niter, rootsol.objval))
    niter += 1
    cut_added = false
    itertime = @mytime rootsol, stats, z_UB, σ = iteration(root, mccount, num_stages, verbose, pathsel, TOL)
    #Lower bound since θ >= 0
    z_LB = rootsol.objval

    totaltime  += itertime
    totalstats += stats
    if verbose >= 2
      println("Iteration $niter completed in $itertime s (Total time is $totaltime)")
      println("Status: $(rootsol.status)")
      println("z_UB: $(z_UB)")
      println("z_LB: $(z_LB)")
      if mccount != :All
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
    if mccount != :All
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
