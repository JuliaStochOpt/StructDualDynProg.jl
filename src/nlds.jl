using MathProgBase

# Primal
# min c x
#     + θ           -> :AveragedCut
#     + sum θ_i     -> :MultiCut
#     Wx = h - Tx_a
#     Dx >= d (feasibility cuts)
#     Ex + θ   >= e -> :AveragedCut
#     Ex + θ_i >= e -> :MultiCut
#     x in K
#     (θ >= 0)

# Dual
# max π (h - Tx_a) + σ d + ρ e
#     c - (π W + σ D + ρ E) in K^*
#     1'ρ = 1 (or <= 1 if θ >= 0)
#     σ >= 0
#     ρ >= 0
type NLDS{S}
  W::AbstractMatrix{S}
  h::Vector{S}
  T::AbstractMatrix{S}
  K
  C
  c::Vector{S}

  # used to generate cuts
  cuts_d
  cuts_e

  # parent solution
  x_a

  childFC::Vector{CutStore{S}}
  childOC::Vector{CutStore{S}}
  localFC::CutStore{S}
  localOC::CutStore{S}
  proba

  nθ
  nπ
  nσ
  nρ
  πs
  σs
  ρs

  model
  loaded

  function NLDS(W::AbstractMatrix{S}, h::Vector{S}, T::AbstractMatrix{S}, K, C, c::Vector{S}, solver)
    nθ = 0
    nπ = length(h)
    localFC = CutStore{S}(size(W, 2))
    localOC = CutStore{S}(size(W, 2))
    if false
      model = MathProgBase.ConicModel(solver)
    else
      model = MathProgBase.LinearQuadraticModel(solver)
    end
    nlds = new(W, h, T, K, C, c, nothing, nothing, S[], CutStore{S}[], CutStore{S}[], localFC, localOC, nothing, nθ, nπ, 0, 0, 1:nπ, nothing, nothing, model, false)
    addfollower(localFC, nlds)
    addfollower(localOC, nlds)
    nlds
  end
end

function setchildren!(nlds::NLDS, childFC, childOC, proba, cutmode)
  @assert length(childFC) == length(childOC) == length(proba)
  if cutmode == :MultiCut
    nlds.proba = proba
    nlds.nθ = length(proba)
  elseif cutmode == :AveragedCut
    nlds.nθ = 1
  else
    nlds.nθ = 0
  end
  nlds.childFC = childFC
  for store in childFC
    addfollower(store, nlds)
  end
  nlds.childOC = childOC
  for store in childOC
    addfollower(store, nlds)
  end
end

function notifynewcut(nlds::NLDS)
  # Invalidate current model
  nlds.loaded = false
end

function setparentx(nlds::NLDS, x_a)
  nlds.x_a = x_a
  # TODO hotstart
  nlds.loaded = false
end

function getfeasibilitycuts(nlds::NLDS)
  cuts_D = reduce(vcat, nlds.localFC.A, map(x -> x.A, nlds.childFC))
  cuts_d = reduce(vcat, nlds.localFC.b, map(x -> x.b, nlds.childFC))
  (cuts_D, cuts_d)
end

function getoptimalitycuts{S}(nlds::NLDS{S})
  function f(i)
    E = nlds.childOC[i].A
    nrows = size(E, 1)
    [E spzeros(S, nrows, i-1) ones(S, nrows, 1) spzeros(S, nrows, nlds.nθ-i)]
  end
  if nlds.nθ == 1
    cuts_E = [nlds.localOC.A ones(size(nlds.localOC.A, 1), 1)]
  else
    cuts_E = reduce(vcat, spzeros(S, 0, size(nlds.W, 2) + nlds.nθ), map(f, 1:length(nlds.childOC)))
  end
  cuts_e = reduce(vcat, nlds.localOC.b, map(x -> x.b, nlds.childOC))
  (cuts_E, cuts_e)
end

function myload!(model::MathProgBase.AbstractConicModel, c, A, b, K, C)
  MathProgBase.loadproblem!(nlds.model, c, A, b, K, C)
end

function myload!(model::MathProgBase.AbstractLinearQuadraticModel, c, A, b, K, C)
  l  = Vector{Float64}(size(A, 2))
  u  = Vector{Float64}(size(A, 2))
  lb = Vector{Float64}(size(A, 1))
  ub = Vector{Float64}(size(A, 1))
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

function load!{S}(nlds::NLDS{S})
  if !nlds.loaded
    bigA = nlds.W
    bigb = nlds.h - nlds.T * nlds.x_a
    nlds.nπ = length(nlds.h)
    cuts_D, nlds.cuts_d = getfeasibilitycuts(nlds)
    bigA = [bigA; cuts_D]
    bigb = [bigb; nlds.cuts_d]
    nlds.nσ = length(nlds.cuts_d)
    if nlds.nθ > 0
      bigA = [bigA spzeros(size(bigA, 1), nlds.nθ)]
      cuts_E, nlds.cuts_e = getoptimalitycuts(nlds)
      bigA = [bigA; cuts_E]
      bigb = [bigb; nlds.cuts_e]
      nlds.nρ = length(nlds.cuts_e)
    else
      nlds.cuts_e = S[]
      nlds.nρ = 0
    end
    nlds.πs = 1:nlds.nπ
    nlds.σs = nlds.nπ+(1:nlds.nσ)
    nlds.ρs = nlds.nπ+nlds.nσ+(1:nlds.nρ)
    if size(bigA, 1) > size(nlds.W, 1)
      bigK = [nlds.K;
             (:NonPos, collect(nlds.σs));
             (:NonPos, collect(nlds.ρs))]
    else
      bigK = nlds.K
    end
    if nlds.nθ == 0
      bigC = nlds.C
      bigc = nlds.c
    else
      bigC = [nlds.C; (:NonNeg, collect(size(nlds.W, 2)+(1:nlds.nθ)))]
      bigc = nlds.c
      if nlds.proba === nothing
        @assert nlds.nθ == 1
        bigc = [bigc; 1]
      else
        bigc = [bigc; nlds.proba]
      end
    end

    myload!(nlds.model, bigc, bigA, bigb, bigK, bigC)
    nlds.loaded = true
  end
end

# π, σ and ρ do not really make sense alone so only
# their product will T, h, d, e is stored
type NLDSSolution
  status
  objval
  x
  θ
  πT
  πh
  σd
  ρe
end

function solve!(nlds::NLDS)
  load!(nlds)
  MathProgBase.optimize!(nlds.model)
  status = MathProgBase.status(nlds.model)
  objval = MathProgBase.getobjval(nlds.model)
  primal = MathProgBase.getsolution(nlds.model)
  dual   = mygetdual(nlds.model)
  x = primal[1:end-nlds.nθ]
  θ = primal[end-nlds.nθ+1:end]
  π = dual[nlds.πs]
  σ = dual[nlds.σs]
  ρ = dual[nlds.ρs]
  NLDSSolution(status, objval, x, θ, vec(π' * nlds.T), vecdot(π, nlds.h), vecdot(σ, nlds.cuts_d), vecdot(ρ, nlds.cuts_e))
end
