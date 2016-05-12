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

function iteration(root::SDDPNode, pathss, num_stages, cutmode=:MultiCut, mccount=25, verbose=0, TOL=1e-5)
  rootsol = nothing

  stats = SDDPStats()
  for t in 1:num_stages
    if verbose >= 3
      @show t
    end
    # children are at t, parents are at t-1
    newpathss = []
    for (parent, psol, paths) in pathss
      if parent == nothing
        rootsol = loadAndSolve(root)
        push!(newpathss, (root, rootsol, paths))
      else
        curpathss = []
        childsolved = zeros(Bool, length(parent.children))
        feasible = true
        nnewfcuts = 0
        nnewocuts = 0
        childocuts = Array{Any}(length(parent.children))
        for i in 1:length(parent.children)
          child = parent.children[i]
          childnpaths = numberofpaths(child, t, num_stages)
          newpaths, paths = filter(childnpaths, paths)
          paths -= childnpaths
          if t == 2 || length(newpaths) > 0 || cutmode == :AveragedCut
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
              feasible = false
              nnewfcuts += 1
              stats.fcutstime += @mytime pushfeasibilitycut!(child, coef, rhs)
              break
            else
              rhs += childsol.ρe
              childocuts[i] = (coef, rhs)
            end
            if length(newpaths) > 0
              push!(curpathss, (child, childsol, newpaths))
            end
          end
        end
        if feasible
          if cutmode == :MultiCut
            for i in 1:length(parent.children)
              if childsolved[i]
                a, β = childocuts[i][1], childocuts[i][2]
                if psol.θ[i] < (β - dot(a, psol.x)) - TOL
                  stats.ocutstime += @mytime pushoptimalitycutforparent!(parent.children[i], a, β)
                  nnewocuts += 1
                end
              end
            end
          elseif cutmode == :AveragedCut
            if !isempty(parent.children)
              # FIXME should use childT
              a = sum(map(i->childocuts[i][1]*parent.proba[i], 1:length(parent.children)))
              β = sum(map(i->childocuts[i][2]*parent.proba[i], 1:length(parent.children)))
              if psol.θ[1] < (β - dot(a, psol.x)) - TOL
                stats.ocutstime += @mytime pushoptimalitycut!(parent, a, β)
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
    end
    pathss = newpathss
  end
  rootsol, stats
end

function SDDP(root::SDDPNode, num_stages, cutmode=:MultiCut, mccount=25, verbose=0, TOL=1e-5)
  # If the graph is not a tree, this will loop if I don't use a num_stages limit
  npaths = numberofpaths(root, 1, num_stages)
  if mccount == :All
    mccount = npaths
  else
    mccount = min(mccount, npaths)
  end

  cut_added = true
  niter = 0
  nfcuts = 0
  nocuts = 0
  if verbose >= 3
    @show npaths
    @show mccount
  end
  rootsol = nothing
  totaltime = 0
  totalstats = SDDPStats()
  while (mccount < npaths || cut_added) && (rootsol === nothing || rootsol.status != :Infeasible)
    niter += 1
    cut_added = false
    if mccount == npaths
      pathss = [(nothing, Float64[], collect(1:npaths))]
    else
      pathss = [(nothing, Float64[], sort(1+mod(rand(Int, mccount), npaths)))]
    end
    itertime = @mytime rootsol, stats = iteration(root, pathss, num_stages, cutmode, mccount, verbose, TOL)
    totaltime  += itertime
    totalstats += stats
    if verbose >= 2
      println("Iteration $niter completed in $itertime s (Total time is $totaltime)")
      println("Status: $(rootsol.status)")
      println("Objective value: $(rootsol.objval)")
      println(" Solution value: $(rootsol.x)")
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
    println("Objective value: $(rootsol.objval)")
    println(" Solution value: $(rootsol.x)")
    println("Total stats:")
    println(totalstats)
  end

  attrs = Dict()
  attrs[:niter] = niter
  attrs[:nfcuts] = nfcuts
  attrs[:nocuts] = nocuts
  SDDPSolution(rootsol.status, rootsol.objval, rootsol.x, attrs)
end
