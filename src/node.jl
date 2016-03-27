import MathProgBase

function loadAndSolveConicProblem(c, A, b, K, C, solver)

    # load conic model
    model = MathProgBase.ConicModel(solver)
    MathProgBase.loadproblem!(model, c, A, b, K, C)

    #println("process id $(myid()) started")

    # solve conic model
    MathProgBase.optimize!(model)
    status = MathProgBase.status(model)

    # return status and dual
    #println("process id $(myid()) status $(status)")
    return status, MathProgBase.getobjval(model), MathProgBase.getsolution(model), MathProgBase.getdual(model)
end

type SDDPNode
  W::Matrix{Float64}
  h::Vector{Float64}
  T::Matrix{Float64}
  K
  C
  c
  solver
  parent
  children::Vector{SDDPNode}
  proba::Vector{Float64}
  # Optimality cuts
  cuts_e::Vector{Float64}
  cuts_E::Matrix{Float64}
  cuts_pe::Vector{Float64}
  cuts_pE::Matrix{Float64}
  # Feasibility Cuts
  cuts_d::Vector{Float64}
  cuts_D::Matrix{Float64}
  root::Bool
  leaf::Bool

  status::Symbol
  objval
  x
  θ
  dual

  function SDDPNode(W, h, T, K, C, c, solver, parent = nothing)
    xsize = size(W, 2)
    cuts_e = Vector{Float64}(0)
    cuts_E = Matrix{Float64}(0, xsize)
    cuts_pe = Vector{Float64}(0)
    cuts_pE = Matrix{Float64}(0, parent !== nothing ? size(parent.W, 2) : 0)
    cuts_d = Vector{Float64}(0)
    cuts_D = Matrix{Float64}(0, xsize)
    root = parent == nothing
    new(W, h, T, K, C, c, solver, parent, SDDPNode[], Float64[], cuts_e, cuts_E, cuts_pe, cuts_pE, cuts_d, cuts_D, root, true, :Undefined, nothing, nothing, nothing, nothing)
  end

end

function setchildren!(node, children, proba)
  @assert length(children) == length(proba)
  node.children = children
  node.proba = proba
  node.leaf = false
end

function pushfeasibilitycut!(node, coef, rhs)
  node.cuts_d = [node.cuts_d; rhs]
  node.cuts_D = [node.cuts_D; coef']
end
function pushoptimalitycut!(node, coef, rhs)
  node.cuts_e = [node.cuts_e; rhs]
  node.cuts_E = [node.cuts_E; coef']
end
function pushoptimalitycutforparent!(node, coef, rhs)
  node.cuts_pe = [node.cuts_pe; rhs]
  node.cuts_pE = [node.cuts_pE; coef']
end

function loadAndSolve(node::SDDPNode, cutmode)
  if node.root
    bigb = node.h
  else
    bigb = node.h - node.T * node.parent.x
  end
  bigA = [node.W; node.cuts_D]
  bigb = [bigb; node.cuts_d]
  if node.leaf
    nθ = 0
  else
    if cutmode == :MultiCut
      nθ = length(node.proba)
    else
      nθ = 1
    end
  end
  if !node.leaf
    bigA = [bigA spzeros(size(bigA, 1), nθ)]
    if cutmode == :MultiCut
      cuts_E = Matrix{Float64}(0, size(bigA, 2))
      cuts_e = Float64[]
      for i in 1:nθ
        nrows = size(node.children[i].cuts_pE, 1)
        cuts_E = [cuts_E; node.children[i].cuts_pE spzeros(nrows, i-1) -ones(nrows, 1) spzeros(nrows, nθ-i)]
        cuts_e = [cuts_e; node.children[i].cuts_pe]
      end
      bigA = [bigA; cuts_E]
      bigb = [bigb; cuts_e]
    else
      bigA = [bigA; node.cuts_E -ones(size(node.cuts_E, 1), 1)]
      bigb = [bigb; node.cuts_e]
    end
  end
  if size(bigA, 1) > size(node.W, 1)
    bigK = [node.K; (:NonNeg, collect((size(node.W, 1)+1):size(bigA, 1)))]
  else
    bigK = node.K
  end
  if node.leaf
    bigC = node.C
    bigc = node.c
  else
    bigC = [node.C; (:NonNeg, collect(size(node.W, 2)+(1:nθ)))]
    if cutmode == :MultiCut
      bigc = [node.c; node.proba]
    else
      bigc = [node.c; 1]
    end
  end
  node.status, node.objval, primal, node.dual = loadAndSolveConicProblem(bigc, bigA, bigb, bigK, bigC, node.solver)
  node.x = primal[1:end-nθ]
  node.θ = primal[end-nθ+1:end]
end

function addCuttingPlanes(node, cutmode, TOL)
  cut_added = false
  infeasible_master = false
  xsize = size(node.W, 2)
  if cutmode == :AveragedCut
    averaged_optimality_cut_coef = zeros(Float64, xsize)
    averaged_optimality_cut_rhs = .0
  end
  # add cutting planes, one per scenario
  for i in 1:length(node.children)
    child = node.children[i]
    coef = vec(child.dual' * child.T)
    rhs = vecdot(child.dual, child.h)
    # add an infeasibility cut
    if child.status == :Infeasible
      # child.dual is a ray
      # so alpha * child.dual is also valid for any alpha >= 0.
      # Hence child.dual might have very large coefficients and alter
      # the numerial accuracy of the master's solver.
      # We scale it to avoid this issue
      scaling = abs(rhs)
      if scaling == 0
        scaling = maximum(coef)
      end
      cut_added = true
      infeasible_master = true
      pushfeasibilitycut!(node, coef/scaling, sign(rhs))
      # add an optimality cut
    else
      if !infeasible_master
        if cutmode == :AveragedCut
          averaged_optimality_cut_coef += node.proba[i] * coef
          averaged_optimality_cut_rhs += node.proba[i] * rhs
        else
          if !infeasible_master && node.θ[i] < (dot(coef, node.x) - rhs) - TOL
            pushoptimalitycutforparent!(child, coef, rhs)
            cut_added = true
          end
        end
      end
    end
  end
  if cutmode == :AveragedCut && !infeasible_master && node.θ[1] < (dot(averaged_optimality_cut_coef, node.x) - averaged_optimality_cut_rhs) - TOL
    pushoptimalitycut!(node, averaged_optimality_cut_coef, averaged_optimality_cut_rhs)
    cut_added = true
  end
  cut_added
end
