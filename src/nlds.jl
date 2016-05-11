export NLDS

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

# Nested L-Shaped Decomposition Subproblem (NLDS)

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
  cuts_de

  # parent solution
  x_a

  childFC::Vector{CutStore{S}}
  childOC::Vector{CutStore{S}}
  localFC::CutStore{S}
  localOC::CutStore{S}
  proba
  childT::Nullable{Vector{AbstractMatrix{S}}}

  nx
  nθ
  nπ
  nσ
  nρ
  πs
  σs
  ρs

  model
  loaded
  solved
  constrsadded
  boundsupdated
  prevsol::Nullable{NLDSSolution}

  newcut::Symbol

  function NLDS(W::AbstractMatrix{S}, h::AbstractVector{S}, T::AbstractMatrix{S}, K, C, c::AbstractVector{S}, solver, newcut::Symbol=:AddImmediately)
    nx = size(W, 2)
    nθ = 0
    nπ = length(h)
    localFC = CutStore{S}(size(W, 2))
    localOC = CutStore{S}(size(W, 2))
    if false
      model = MathProgBase.ConicModel(solver)
    else
      model = MathProgBase.LinearQuadraticModel(solver)
    end
    nlds = new(W, h, T, K, C, c, nothing, S[], CutStore{S}[], CutStore{S}[], localFC, localOC, nothing, nothing, nx, nθ, nπ, 0, 0, 1:nπ, nothing, nothing, model, false, false, false, false, nothing, newcut)
    addfollower(localFC, (nlds, (:Feasibility, 0)))
    addfollower(localOC, (nlds, (:Optimality, 0)))
    nlds
  end
end

function (::Type{NLDS{S}}){S}(W::AbstractMatrix, h::AbstractVector, T::AbstractMatrix, K, C, c::AbstractVector, solver, newcut::Symbol=:AddImmediately)
  NLDS{S}(AbstractMatrix{S}(W), AbstractVector{S}(h), AbstractMatrix{S}(T), K, C, AbstractVector{S}(c), solver, newcut)
end

function setchildren!(nlds::NLDS, childFC, childOC, proba, cutmode, childT=nothing)
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
  for i in 1:length(childFC)
    addfollower(childFC[i], (nlds, (:Feasibility, i)))
  end
  nlds.childOC = childOC
  for i in 1:length(childOC)
    addfollower(childOC[i], (nlds, (:Optimality, i)))
  end
  nlds.childT = childT
end

function getfeasibilitycuts(nlds::NLDS)
  function f(i)
    D = nlds.childFC[i].A
    if !isnull(nlds.childT)
      D = D * get(nlds.childT)[i]
    end
    D
  end
  cuts_D = reduce(vcat, nlds.localFC.A, map(i -> f(i), 1:length(nlds.childFC)))
  cuts_d = reduce(vcat, nlds.localFC.b, map(x -> x.b, nlds.childFC))
  (cuts_D, cuts_d)
end

function getoptimalitycuts{S}(nlds::NLDS{S})
  function f(i)
    E = nlds.childOC[i].A
    if !isnull(nlds.childT)
      E = E * get(nlds.childT)[i]
    end
    nrows = size(E, 1)
    [E spzeros(S, nrows, i-1) ones(S, nrows, 1) spzeros(S, nrows, nlds.nθ-i)]
  end
  if nlds.nθ == 1
    cuts_E = [nlds.localOC.A ones(size(nlds.localOC.A, 1), 1)]
  else
    cuts_E = spzeros(S, 0, nlds.nx + nlds.nθ)
  end
  if nlds.nθ == length(nlds.childOC)
    cuts_E = reduce(vcat, spzeros(S, 0, nlds.nx + nlds.nθ), map(f, 1:length(nlds.childOC)))
  end
  cuts_e = reduce(vcat, nlds.localOC.b, map(x -> x.b, nlds.childOC))
  (cuts_E, cuts_e)
end

function notifynewcut{S}(nlds::NLDS{S}, a::AbstractVector{S}, β::S, attrs)
  if nlds.loaded
    if nlds.newcut == :InvalidateSolver
      nlds.loaded = false
      nlds.solved = false
    elseif nlds.newcut == :AddImmediately
      @assert attrs[1] in [:Feasibility, :Optimality]
      idx = collect(1:nlds.nx)
      if attrs[1] == :Feasibility
        nlds.nσ += 1
        push!(nlds.σs, nlds.nπ + nlds.nσ + nlds.nρ)
      else
        a = [a; one(S)]
        if attrs[2] == 0
          push!(idx, nlds.nx+1)
        else
          push!(idx, nlds.nx+attrs[2])
        end
        nlds.nρ += 1
        push!(nlds.ρs, nlds.nπ + nlds.nσ + nlds.nρ)
      end
      if attrs[2] > 0 && !isnull(nlds.childT)
        a = get(nlds.childT)[attrs[2]]' * a
      end
      applyboundsupdates!(nlds)
      myaddconstr!(nlds.model, idx, a, β, :NonPos)
      #push!(nlds.cuts_de, β)
      nlds.cuts_de = [nlds.cuts_de; sparsevec([β])]
      nlds.constrsadded = true
      nlds.solved = false
    else
      error("Invalid newcut option $(nlds.newcut)")
    end
  end
end

function getrhs(nlds)
  bigb = [nlds.h - nlds.T * nlds.x_a; nlds.cuts_de]
  bigK = nlds.K
  if nlds.nσ > 0
    push!(bigK, (:NonPos, nlds.σs))
  end
  if nlds.nρ > 0
    push!(bigK, (:NonPos, nlds.ρs))
  end
  bigb, bigK
end

function setparentx(nlds::NLDS, x_a)
  nlds.x_a = x_a
  if nlds.loaded
    bigb, bigK = getrhs(nlds)
    applyconstraintsaddition!(nlds)
    mysetconstrB!(nlds.model, bigb, bigK)
    nlds.boundsupdated = true
    nlds.solved = false
  end
end

function applyboundsupdates!(nlds::NLDS)
  if nlds.boundsupdated
    MathProgBase.updatemodel!(nlds.model)
    nlds.boundsupdated = false
  end
end

function applyconstraintsaddition!(nlds::NLDS)
  if nlds.constrsadded
    MathProgBase.updatemodel!(nlds.model)
    nlds.constrsadded = false
  end
end

function load!{S}(nlds::NLDS{S})
  if !nlds.loaded
    bigA = nlds.W
    nlds.nπ = length(nlds.h)
    cuts_D, cuts_d = getfeasibilitycuts(nlds)
    bigA = [bigA; cuts_D]
    nlds.nσ = length(cuts_d)
    if nlds.nθ > 0
      bigA = [bigA spzeros(size(bigA, 1), nlds.nθ)]
      cuts_E, cuts_e = getoptimalitycuts(nlds)
      bigA = [bigA; cuts_E]
      nlds.nρ = length(cuts_e)
    else
      cuts_e = S[]
      nlds.nρ = 0
    end
    nlds.cuts_de = [cuts_d; cuts_e]
    nlds.πs = collect(1:nlds.nπ)
    nlds.σs = collect(nlds.nπ+(1:nlds.nσ))
    nlds.ρs = collect(nlds.nπ+nlds.nσ+(1:nlds.nρ))
    if nlds.nθ == 0
      bigC = nlds.C
      bigc = nlds.c
    else
      bigC = [nlds.C; (:NonNeg, collect(nlds.nx+(1:nlds.nθ)))]
      bigc = nlds.c
      if nlds.proba === nothing
        @assert nlds.nθ == 1
        bigc = [bigc; 1]
      else
        bigc = [bigc; nlds.proba]
      end
    end

    bigb, bigK = getrhs(nlds)

    myload!(nlds.model, bigc, bigA, bigb, bigK, bigC)
    nlds.loaded = true
  end
end

function solve!(nlds::NLDS)
  load!(nlds)
  if !nlds.solved
    applyboundsupdates!(nlds)
    # No need to call updatemodel!(model) for added constraints for Gurobi and only Gurobi needs model updates
    nlds.constrsadded = false
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
    cuts_d = nlds.cuts_de[nlds.σs-nlds.nπ]
    cuts_e = nlds.cuts_de[nlds.ρs-nlds.nπ]
    nlds.prevsol = NLDSSolution(status, objval, x, θ, vec(π' * nlds.T), vecdot(π, nlds.h), vecdot(σ, cuts_d), vecdot(ρ, cuts_e))
  end
end

function getsolution(nlds::NLDS)
  solve!(nlds)
  get(nlds.prevsol)
end
