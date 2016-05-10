type SDDPSolution
    status
    objval
    sol
    attrs
end

function SDDP(root::SDDPNode, num_stages, cutmode=:MultiCut, mccount=25, debug=false, TOL=1e-5)
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
  if debug
    @show npaths
    @show mccount
  end
  rootsol = nothing
  while (mccount < npaths || cut_added) && (rootsol === nothing || rootsol.status != :Infeasible)
    if debug
      @show niter
    end
    niter += 1
    cut_added = false
    pathss = [(nothing, Float64[], sort(randperm(npaths)[1:mccount]))]
    for t in 1:num_stages
      if debug
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
            if length(newpaths) > 0 || cutmode == :AveragedCut
              setchildx(parent, i, psol.x)
              childsol = loadAndSolve(child)
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
                pushfeasibilitycut!(child, coef, rhs)
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
                    pushoptimalitycutforparent!(parent.children[i], a, β)
                    nnewocuts += 1
                  end
                end
              end
            elseif cutmode == :AveragedCut
              if !isempty(parent.children)
                a = sum(map(i->childocuts[i][1]*parent.proba[i], 1:length(parent.children)))
                β = sum(map(i->childocuts[i][2]*parent.proba[i], 1:length(parent.children)))
                if psol.θ[1] < (β - dot(a, psol.x)) - TOL
                  pushoptimalitycut!(parent, a, β)
                  nnewocuts += 1
                end
              end
            end
          end
          nfcuts += nnewfcuts
          nocuts += nnewocuts
          cut_added |= (nnewfcuts + nnewocuts > 0)
          if feasible && !isempty(curpathss)
            @assert nnewfcuts == 0
            append!(newpathss, curpathss)
          else
            @assert nnewfcuts == 1
          end
        end
      end
      pathss = newpathss
    end
    if debug
      @show rootsol.status, rootsol.objval, rootsol.x, niter, nfcuts, nocuts
    end
  end

  attrs = Dict()
  attrs[:niter] = niter
  attrs[:nfcuts] = nfcuts
  attrs[:nocuts] = nocuts
  SDDPSolution(rootsol.status, rootsol.objval, rootsol.x, attrs)
end
