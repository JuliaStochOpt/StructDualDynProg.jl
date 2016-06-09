export waitandsee
import Gurobi

function waitandsee(root::SDDPNode, num_stages, solver, totalmccount=25, verbose=0)
  typealias Path Tuple{SDDPNode, Vector{NLDS}, Float64, Float64, typeof(totalmccount)}
  paths = Path[(root, NLDS[root.nlds], .0, 1., totalmccount)]
  for t in 2:num_stages
    newpaths = Path[]
    for (parent, nldspath, z, prob, mccount) in paths
      if isempty(parent.children)
        push!(newpaths, (parent, nldspath, z, prob, mccount))
      else
        npaths = choosepaths(parent, mccount, :Proba, t, num_stages)
        childs = totalmccount == :All ? (1:length(parent.children)) : find(npaths .> 0)
        curpaths = Path[(parent.children[i], [nldspath; parent.children[i].nlds], z, prob*parent.proba[i], npaths[i]) for i in childs]
        append!(newpaths, curpaths)
      end
    end
    paths = newpaths
  end

  sumz = 0
  newpaths = Path[]
  for (last, nldspath, z, prob, mccount) in paths
    nvars = cumsum(Int[nlds.nx for nlds in nldspath])
    ncons = cumsum(Int[nlds.nπ for nlds in nldspath])
    # b - Ax in K_1
    # x \in K_2
    A = spzeros(ncons[end], nvars[end])
    bs = Vector{Vector}(length(nldspath))
    Ks = []
    C = []
    c = similar(root.nlds.c, nvars[end])
    for i in 1:length(nldspath)
      offsetnvars = i == 1 ? 0 : nvars[i-1]
      offsetncons = i == 1 ? 0 : ncons[i-1]
      nlds = nldspath[i]
      c[offsetnvars+1:nvars[i]] = nlds.c
      bs[i] = nlds.h
      A[offsetncons+1:ncons[i], offsetnvars+1:nvars[i]] = nlds.W
      if i > 1
        offset = i == 2 ? 0 : nvars[i-2]
        A[offsetncons+1:ncons[i], offset+1:nvars[i-1]] = nlds.T
      end
      push!(Ks, nlds.K)
      append!(C, [(cone, offsetnvars+idx) for (cone, idx) in nlds.C])
    end
    model = MathProgBase.LinearQuadraticModel(solver)
    myload!(model, c, A, bs, Ks, C)
    MathProgBase.optimize!(model)
    status = MathProgBase.status(model)
    if status == :Error
      error("The solver reported an error")
    elseif status == :UserLimit
      error("The solver reached iteration limit or timed out")
    elseif status == :Infeasible
      error("The problem is infeasible for some realization of the uncertainty")
    elseif status == :Unbounded
      error("The problem is unbounded for some realization of the uncertainty")
    else
      @assert status == :Optimal
      objval = MathProgBase.getobjval(model)
      push!(newpaths, (last, nldspath, objval, prob, mccount))
    end
  end
  meanstdpaths(newpaths, 3, 4, 5, totalmccount)
end
