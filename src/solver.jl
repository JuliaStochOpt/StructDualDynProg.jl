using MathProgBase

function mysetconstrB!(m::MathProgBase.AbstractConicModel, b, K)
  error("Not supported")
end

function getLPconstrbounds(b, K)
  lb = Vector{Float64}(length(b))
  ub = Vector{Float64}(length(b))
  for (cone, idx) in K
    for i in idx
      if !(cone  in [:Zero, :NonPos, :NonNeg])
        error("This cone is not supported")
      end
      if cone == :Zero || cone == :NonPos
        lb[i] = b[i]
      else
        lb[i] = -Inf
      end
      if cone == :Zero || cone == :NonNeg
        ub[i] = b[i]
      else
        ub[i] = Inf
      end
    end
  end
  lb, ub
end

function mysetconstrB!(m::MathProgBase.AbstractLinearQuadraticModel, b, K)
  lb, ub = getLPconstrbounds(b, K)
  MathProgBase.setconstrLB!(m, lb)
  MathProgBase.setconstrUB!(m, ub)
end

function myaddconstr!(m::MathProgBase.AbstractConicModel, idx, a, β, cone)
  error("Not supported")
end

function myaddconstr!(m::MathProgBase.AbstractLinearQuadraticModel, idx, a, β, cone)
  lb = -Inf
  ub = Inf
  if cone == :Zero || cone == :NonPos
    lb = β
  end
  if cone == :Zero || cone == :NonNeg
    ub = β
  end
  MathProgBase.addconstr!(m, idx, a, lb, ub)
end

function myload!(model::MathProgBase.AbstractConicModel, c, A, b, K, C)
  MathProgBase.loadproblem!(nlds.model, c, A, b, K, C)
end

function myload!(model::MathProgBase.AbstractLinearQuadraticModel, c, A, b, K, C)
  l  = Vector{Float64}(size(A, 2))
  u  = Vector{Float64}(size(A, 2))
  lb, ub = getLPconstrbounds(b, K)
  for (cone, idx) in C
    for i in idx
      if !(cone  in [:Free, :NonPos, :NonNeg])
        error("This cone is not supported")
      end
      if cone == :Free || cone == :NonPos
        l[i] = -Inf
      else
        l[i] = 0
      end
      if cone == :Free || cone == :NonNeg
        u[i] = Inf
      else
        u[i] = 0
      end
    end
  end
  MathProgBase.loadproblem!(model, A, l, u, c, lb, ub, :Min)
end

function mygetdual(model::MathProgBase.AbstractConicModel)
  -MathProgBase.getdual(model)
end
function mygetdual(model::MathProgBase.AbstractLinearQuadraticModel)
  if MathProgBase.status(model) == :Infeasible
    MathProgBase.getinfeasibilityray(model)
  else
    MathProgBase.getconstrduals(model)
  end
end
