import MathProgBase

export SDDPNode, setchildren!

type SDDPNode{S}
  nlds::NLDS{S}
  nvars::Int
  parent::Nullable{SDDPNode{S}}
  children::Vector{SDDPNode{S}}
  proba::Vector{S}
  # Optimality cuts
  root::Bool
  leaf::Bool
  npath::Nullable{Int}

  fcuts::CutStore{S}
  ocuts::CutStore{S}

  sol

  function SDDPNode(nlds::NLDS{S}, parent)
    nvars = size(nlds.W, 2)
    root = parent === nothing
    nvars_a = root ? 0 : parent.nvars
    new(nlds, nvars, parent, SDDPNode[], Float64[], root, true, nothing, CutStore{S}(nvars_a), CutStore{S}(nvars_a), nothing)
  end

end

SDDPNode{S}(nlds::NLDS{S}, parent) = SDDPNode{S}(nlds, parent)

function setchildren!(node::SDDPNode, children, proba, cutmode)
  @assert length(children) == length(proba)
  node.children = children
  node.proba = proba
  node.leaf = false
  childFC = map(child -> child.fcuts, children)
  childOC = map(child -> child.ocuts, children)
  setchildren!(node.nlds, childFC, childOC, proba, cutmode)
end

function numberofpaths(node::SDDPNode)
  if isnull(node.npath)
    if length(node.children) == 0
      node.npath = 1
    else
      node.npath = sum(map(numberofpaths, node.children))
    end
  end
  get(node.npath)
end

function pushfeasibilitycut!(node, coef, rhs)
  addcut(node.fcuts, coef, rhs)
end
function pushoptimalitycut!(node, coef, rhs)
  addcut(node.nlds.localOC, coef, rhs)
end
function pushoptimalitycutforparent!(node, coef, rhs)
  addcut(node.ocuts, coef, rhs)
end

# Feasibility cut
# D = π T
# d = π h + σ d
# Optimality cut
# E = π T
# e = π h + ρ e + σ d
function addCuttingPlanes(node, cutmode, TOL)
  nfcuts = 0
  nocuts = 0
  infeasible_master = false
  if cutmode == :AveragedCut
    averaged_optimality_cut_coef = zeros(Float64, node.nvars)
    averaged_optimality_cut_rhs = .0
  end
  # add cutting planes, one per scenario
  for i in 1:length(node.children)
    child = node.children[i]
    coef = child.sol.πT
    rhs = child.sol.πh + child.sol.σd
    # add an infeasibility cut
    if child.sol.status == :Infeasible
      # child.dual is a ray
      # so alpha * child.dual is also valid for any alpha >= 0.
      # Hence child.dual might have very large coefficients and alter
      # the numerial accuracy of the master's solver.
      # We scale it to avoid this issue
      scaling = abs(rhs)
      if scaling == 0
        scaling = maximum(coef)
      end
      infeasible_master = true
      pushfeasibilitycut!(child, coef/scaling, sign(rhs))
      #pushfeasibilitycut!(node, coef, rhs)
      nfcuts += 1
      # add an optimality cut
    else
      rhs += child.sol.ρe
      if !infeasible_master
        if cutmode == :AveragedCut
          averaged_optimality_cut_coef += node.proba[i] * coef
          averaged_optimality_cut_rhs += node.proba[i] * rhs
        else
          if !infeasible_master && node.sol.θ[i] < (rhs - dot(coef, node.sol.x)) - TOL
            pushoptimalitycutforparent!(child, coef, rhs)
            # FIXME it will be added multiple times ! each time a parent satisfies theta < ...
            nocuts += 1
          end
        end
      end
    end
  end
  if cutmode == :AveragedCut && !infeasible_master && node.sol.θ[1] < (averaged_optimality_cut_rhs - dot(averaged_optimality_cut_coef, node.sol.x)) - TOL
    pushoptimalitycut!(node, averaged_optimality_cut_coef, averaged_optimality_cut_rhs)
    nocuts += 1
  end
  (nfcuts, nocuts)
end

function loadAndSolve(node::SDDPNode)
# if !isnull(node.parent)
#   setparentx(node.nlds, get(node.parent).sol.x)
# end
  node.sol = solve!(node.nlds)
end
