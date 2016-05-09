type SDDPSolution
    status
    objval
    sol
    attrs
end

function SDDP(root::SDDPNode, num_stages, cutmode=:MultiCut, mccount=25, TOL=1e-5)
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
  @show npaths
  @show mccount
  while (mccount < npaths || cut_added) && (root.sol === nothing || root.sol.status != :Infeasible)
    @show niter
    niter += 1
    cut_added = false
    pathss = [(nothing, Float64[], sort(randperm(npaths)[1:mccount]))]
    for t in 1:num_stages
      @show t
      # children are at t, parents are at t-1
      newpathss = []
      for (parent, x, paths) in pathss
        if parent == nothing
          loadAndSolve(root)
          push!(newpathss, (root, root.sol.x, paths))
        else
          curpathss = []
          for child in parent.children
            setparentx(child.nlds, x)
            loadAndSolve(child)
            childnpaths = numberofpaths(child, t, num_stages)
            newpaths, paths = filter(childnpaths, paths)
            paths -= childnpaths
            if length(newpaths) > 0
              push!(curpathss, (child, child.sol.x, newpaths))
            end
          end
          (nnewfcuts, nnewocuts) = addCuttingPlanes(parent, cutmode, TOL)
          nfcuts += nnewfcuts
          nocuts += nnewocuts
          cut_added |= (nnewfcuts + nnewocuts > 0)
          if nnewfcuts == 0 && !isempty(curpathss)
            # Do not continue of there were feasibility cuts
            append!(newpathss, curpathss)
          end
        end
      end
      pathss = newpathss
    end
    @show root.sol.status, root.sol.objval, root.sol.x, niter, nfcuts, nocuts
  end

  attrs = Dict()
  attrs[:niter] = niter
  attrs[:nfcuts] = nfcuts
  attrs[:nocuts] = nocuts
  SDDPSolution(root.sol.status, root.sol.objval, root.sol.x, attrs)
end
