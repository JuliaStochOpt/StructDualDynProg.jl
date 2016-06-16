export waitandsee

type WSPath
  node::SDDPNode
  nlds::Vector{NLDS}
  z::Float64
  proba::Float64
  mccount::Int
end

function waitandsee(root::SDDPNode, num_stages, solver, totalmccount=25, verbose=0)
  paths = WSPath[WSPath(root, NLDS[root.nlds], .0, 1., totalmccount)]
  for t in 2:num_stages
    newpaths = WSPath[]
    for path in paths
      if isempty(path.node.children)
        push!(newpaths, path)
      else
        npaths = choosepaths(path.node, path.mccount, :Proba, t, num_stages)
        childs = totalmccount == -1 ? (1:length(path.node.children)) : find(npaths .> 0)
        curpaths = WSPath[WSPath(path.node.children[i], [path.nlds; path.node.children[i].nlds], path.z, path.proba * path.node.proba[i], npaths[i]) for i in childs]
        append!(newpaths, curpaths)
      end
    end
    paths = newpaths
  end

  sumz = 0
  newpaths = WSPath[]
  for path in paths
    nvars = cumsum(Int[nlds.nx for nlds in path.nlds])
    ncons = cumsum(Int[nlds.nπ for nlds in path.nlds])
    # b - Ax in K_1
    # x \in K_2
    A = spzeros(ncons[end], nvars[end])
    bs = Vector{Vector}(length(path.nlds))
    Ks = []
    C = []
    c = similar(root.nlds.c, nvars[end])
    for i in 1:length(path.nlds)
      offsetnvars = i == 1 ? 0 : nvars[i-1]
      offsetncons = i == 1 ? 0 : ncons[i-1]
      nlds = path.nlds[i]
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
      push!(newpaths, WSPath(path.node, path.nlds, objval, path.proba, path.mccount))
    end
  end
  meanstdpaths(newpaths, totalmccount)
end
