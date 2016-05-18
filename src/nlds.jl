export NLDS, updatemaxncuts!

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

function defaultbettercut(nwitha, nuseda, mycuta, nwithb, nusedb, mycutb)
  # If they are equal, then since we are asked whether "a" is striclty better than "b", we return false.
  # The stability of the sorting algorithm will do the rest.
  # For this reason, only strict inequalities are used ">", "<"

  # I use floating point arithmetic to avoid Int64 overflow
  if mycuta == mycutb
    if nwitha == 0 && nwithb == 0
      false
    elseif nwitha > 0 && nwithb == 0
      if mycutb
        false
      else
        nuseda / nwitha > 3/4
      end
    elseif nwitha == 0 && nwithb > 0
      if mycuta
        true
      else
        nusedb / nwithb < 3/4
      end
    else # nwitha > 0 && nwithb > 0
      nuseda / nwitha > nusedb / nwithb
    end
  elseif mycuta && !mycutb
    if nwitha == 0
      true
    elseif nwithb == 0 # nwitha > 0
      # if the ratio of "a" is larger than 1//2,
      # then we want to keep it over a new cut not created by this NLDS
      nuseda / nwitha > 1/2
    else # nwitha > 0 && nwithb > 0
      # "a" gets a 1/4 bonus trust because it was created by this NLDS
      nuseda / nwitha + 1/4 > nusedb / nwithb
    end
  else
    @assert mycutb && !mycuta
    !defaultbettercut(nwithb, nusedb, mycutb, nwitha, nuseda, mycuta)
  end
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
  cuts_DE::Nullable{AbstractMatrix{S}}
  cuts_de::Nullable{AbstractVector{S}}

  # parent solution
  x_a

  childFC::Vector{CutStore{S}}
  childOC::Vector{CutStore{S}}
  localFC::CutStore{S}
  localOC::CutStore{S}
  proba
  childT::Nullable{Vector{AbstractMatrix{S}}}
  cutmode::Symbol

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
  prevsol::Nullable{NLDSSolution}

  newcut::Symbol
  maxncuts::Integer

  nwith::Vector{Int}
  nused::Vector{Int}
  mycut::Vector{Bool}
  trust::Nullable{Vector{Float64}}

  newcuttrust::Float64
  mycutbonus::Float64
  bettercut::Nullable{Function}

  function NLDS(W::AbstractMatrix{S}, h::AbstractVector{S}, T::AbstractMatrix{S}, K, C, c::AbstractVector{S}, solver, newcut::Symbol=:AddImmediately, maxncuts::Integer=-1, newcuttrust=3/4, mycutbonus=1/4, bettercut=nothing)
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
    nlds = new(W, h, T, K, C, c, nothing, nothing, S[], CutStore{S}[], CutStore{S}[], localFC, localOC, nothing, nothing, :NoOptimalityCut, nx, nθ, nπ, 0, 0, 1:nπ, Int[], Int[], model, false, false, nothing, newcut, maxncuts, Int[], Int[], Bool[], nothing, newcuttrust, mycutbonus, bettercut)
    addfollower(localFC, (nlds, (:Feasibility, 0)))
    addfollower(localOC, (nlds, (:Optimality, 0)))
    nlds
  end
end

function (::Type{NLDS{S}}){S}(W::AbstractMatrix, h::AbstractVector, T::AbstractMatrix, K, C, c::AbstractVector, solver, newcut::Symbol=:AddImmediately, maxncuts::Integer=-1, newcuttrust=3/4, mycutbonus=1/4, bettercut=nothing)
  NLDS{S}(AbstractMatrix{S}(W), AbstractVector{S}(h), AbstractMatrix{S}(T), K, C, AbstractVector{S}(c), solver, newcut, maxncuts, newcuttrust, mycutbonus, bettercut)
end

function setchildren!(nlds::NLDS, childFC, childOC, proba, cutmode, childT)
  nlds.cutmode = cutmode
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

function appendchildren!(nlds::NLDS, childFC, childOC, proba, childT)
  @assert length(childFC) == length(childOC)
  @assert length(proba) == length(nlds.childOC) + length(childOC) == length(nlds.childFC) + length(childFC)
  if nlds.cutmode == :MultiCut
    nlds.proba = proba
    nlds.nθ = length(proba)
  end
  for i in 1:length(childFC)
    @assert length(childFC[i].b) == 0
    addfollower(childFC[i], (nlds, (:Feasibility, length(nlds.childFC)+i)))
  end
  append!(nlds.childFC, childFC)
  for i in 1:length(childOC)
    @assert length(childFC[i].b) == 0
    addfollower(childOC[i], (nlds, (:Optimality, length(nlds.childOC)+i)))
  end
  append!(nlds.childOC, childOC)
  if childT === nothing
    @assert isnull(nlds.childT)
  else
    # If there isn't any child yet, nlds.childT is null
    if isnull(nlds.childT)
      nlds.childT = childT
    else
      append!(get(nlds.childT), childT)
    end
  end
end

function updatemaxncuts!(nlds::NLDS, maxncuts)
  nlds.maxncuts = maxncuts
end

# .=== doesn't work :(
function veceqeqeq(v::Vector, x)
  map(el->el === x, v)
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
  mycut = reduce(vcat, veceqeqeq(nlds.localFC.authors, nlds), map(x -> veceqeqeq(x.authors, nlds), nlds.childFC))
  (cuts_D, cuts_d, mycut)
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
  mycut = reduce(vcat, veceqeqeq(nlds.localOC.authors, nlds), map(x -> veceqeqeq(x.authors, nlds), nlds.childOC))
  (cuts_E, cuts_e, mycut)
end

function gettrustof(nlds::NLDS, nwith, nused, mycut)
  (nwith == 0 ? nlds.newcuttrust : nused / nwith) + (mycut ? nlds.mycutbonus : 0)
end
function updatetrust(nlds::NLDS, i::Int)
  get(nlds.trust)[i] = gettrustof(nlds, nlds.nwith[i], nlds.nused[i], nlds.mycut[i])
end
function gettrust(nlds::NLDS)
  if isnull(nlds.trust)
    trust = nlds.nused ./ nlds.nwith
    trust[nlds.nwith .== 0] = nlds.newcuttrust
    trust[nlds.mycut] += nlds.mycutbonus
    nlds.trust = trust
  end
  get(nlds.trust)
end


function choosecutstoremove(nlds::NLDS, num)
  # MergeSort is stable so in case of equality, the oldest cut loose
  # However PartialQuickSort is a lot faster

  trust = gettrust(nlds)
  if num == 1
    [indmin(trust)]                   # indmin selects the oldest cut in case of tie -> good :)
  else
    sortperm(trust, alg=PartialQuickSort(num))[1:num] # PartialQuickSort is unstable ->  bad :(
  end

# idx = collect(1:length(nlds.nwith))
# sort!(idx, alg=PartialQuickSort(num), lt=(i,j)->nlds.bettercut(nlds.nwith[i], nlds.nused[i], nlds.mycut[i], nlds.nwith[j], nlds.nused[j], nlds.mycut[j]), rev=true)
# idx[1:num]

# isold = nlds.nwith .>= 1
# if reduce(|, false, isold) # at least one old
#   all = 1:length(nlds.nwith)
#   idx = all[isold]
#   if length(idx) > num
#     sorted = sort(idx, by=i->nlds.nused[i]/nlds.nwith[i])
#     sorted[1:num]
#   else
#     new = all[!isold][1:(num-length(idx))]
#     [idx; new]
#   end
# else
#   # No cut has already been used, remove the oldest
#   collect(1:num)
# end
end

function notifynewcut{S}(nlds::NLDS{S}, a::AbstractVector{S}, β::S, attrs, author)
  @assert attrs[1] in [:Feasibility, :Optimality]
  if !isnull(nlds.cuts_DE)
    @assert !isnull(nlds.cuts_de)
    @assert length(nlds.nwith) == length(nlds.nused) == length(nlds.mycut) == length(get(nlds.cuts_de))
    i = attrs[2]
    if i > 0 && !isnull(nlds.childT)
      a = get(nlds.childT)[i]' * a
    end
    if attrs[1] == :Feasibility
      A = [a; spzeros(S, nlds.nθ)]
    else
      if i == 0
        A = [a; one(S)]
      else
        A = [a; spzeros(S, i-1); one(S); spzeros(S, nlds.nθ-i)]
      end
    end
    mine = author === nlds
    if nlds.maxncuts != -1 && length(nlds.nwith) >= nlds.maxncuts
      # Need to remove some cuts
      J = choosecutstoremove(nlds, length(nlds.nwith) - nlds.maxncuts + 1)
      #if nlds.trustnlds.bettercut(0, 0, mine, nlds.nwith[J[end]], nlds.nused[J[end]], nlds.mycut[J[end]])
      if gettrustof(nlds, 0, 0, mine) >= gettrust(nlds)[J[end]]
        j = J[end]
        get(nlds.cuts_DE)[j,:] = sparse(A)
        get(nlds.cuts_de)[j] = β
        nlds.nwith[j] = 0
        nlds.nused[j] = 0
        nlds.mycut[j] = mine
        gettrust(nlds)[j] = gettrustof(nlds, 0, 0, mine)
        cutadded = true
      else
        cutadded = false
      end
      J = J[1:end-1]

      if length(J) > 1 || cutadded
        keep = setdiff(1:length(get(nlds.cuts_de)), J)
        isσcut = zeros(Bool, length(nlds.nwith))
        isσcut[nlds.σs-nlds.nπ] = true
        if cutadded
          isσcut[j] = attrs[1] == :Feasibility
        end
        isσcut = isσcut[keep]
        nlds.σs = nlds.nπ + (1:length(isσcut))[isσcut]
        nlds.ρs = nlds.nπ + (1:length(isσcut))[!isσcut]
        nlds.nσ = length(nlds.σs)
        nlds.nρ = length(nlds.ρs)
      end

      if length(J) > 1
        nlds.cuts_DE = get(nlds.cuts_DE)[keep,:]
        nlds.cuts_de = get(nlds.cuts_de)[keep]
        nlds.nwith = nlds.nwith[keep]
        nlds.nused = nlds.nused[keep]
        nlds.mycut = nlds.mycut[keep]
        nlds.trust = gettrust(nlds)[keep]
      end
      cutremoved = true
    else
      # Just append cut
      if attrs[1] == :Feasibility
        nlds.nσ += 1
        push!(nlds.σs, nlds.nπ + nlds.nσ + nlds.nρ)
      else
        nlds.nρ += 1
        push!(nlds.ρs, nlds.nπ + nlds.nσ + nlds.nρ)
      end
      nlds.cuts_DE = [get(nlds.cuts_DE); sparse(A')]
      nlds.cuts_de = [get(nlds.cuts_de); sparsevec([β])]
      push!(nlds.nwith, 0)
      push!(nlds.nused, 0)
      push!(nlds.mycut, mine)
      if !isnull(nlds.trust)
        push!(get(nlds.trust), gettrustof(nlds, 0, 0, mine))
      end
      cutadded = true
      cutremoved = false
    end
    if nlds.loaded
      if cutremoved || nlds.newcut == :InvalidateSolver
        nlds.loaded = false
        nlds.solved = false
      elseif nlds.newcut == :AddImmediately
        idx = collect(1:nlds.nx)
        if attrs[1] == :Optimality
          if i == 0
            push!(idx, nlds.nx+1)
          else
            push!(idx, nlds.nx+i)
          end
          a = [a; one(S)]
        end
        myaddconstr!(nlds.model, idx, a, β, :NonPos)
        #push!(nlds.cuts_de, β)
        nlds.solved = false
      else
        error("Invalid newcut option $(nlds.newcut)")
      end
    end
  end
end

function checkconsistency(nlds)
  @assert length(nlds.nwith) == length(nlds.nused) == length(nlds.mycut)
  @assert length(nlds.πs) == nlds.nπ
  @assert length(nlds.σs) == nlds.nσ
  @assert length(nlds.ρs) == nlds.nρ
  @assert sort([nlds.πs; nlds.σs; nlds.ρs]) == collect(1:(nlds.nπ + nlds.nσ + nlds.nρ))
end

function getrhs(nlds)
  bigb = [nlds.h - nlds.T * nlds.x_a; get(nlds.cuts_de)]
  # I do a copy since I am going to push!
  bigK = copy(nlds.K)
  if nlds.nσ > 0
    push!(bigK, (:NonPos, nlds.σs))
  end
  if nlds.nρ > 0
    push!(bigK, (:NonPos, nlds.ρs))
  end
  @assert length(bigb) == nlds.nπ + nlds.nσ + nlds.nρ
  bigb, bigK
end

function setparentx(nlds::NLDS, x_a)
  nlds.x_a = x_a
  if nlds.loaded
    bigb, bigK = getrhs(nlds)
    mysetconstrB!(nlds.model, bigb, bigK)
    nlds.solved = false
  end
end

function computecuts!{S}(nlds::NLDS{S})
  if isnull(nlds.cuts_DE)
    @assert isnull(nlds.cuts_de)
    nlds.nπ = length(nlds.h)
    cuts_D, cuts_d, mycut_d = getfeasibilitycuts(nlds)
    nlds.nσ = length(cuts_d)
    if nlds.nθ > 0
      cuts_D = [cuts_D spzeros(S, length(cuts_d), nlds.nθ)]
    end
    #if nlds.nθ > 0
      cuts_E, cuts_e, mycut_e = getoptimalitycuts(nlds)
      nlds.nρ = length(cuts_e)
    #else
    #  cuts_e = S[]
    #  nlds.nρ = 0
    #end
    nlds.cuts_DE = [cuts_D; cuts_E]
    nlds.cuts_de = [cuts_d; cuts_e]
    nlds.πs = collect(1:nlds.nπ)
    nlds.σs = collect(nlds.nπ+(1:nlds.nσ))
    nlds.ρs = collect(nlds.nπ+nlds.nσ+(1:nlds.nρ))
    nlds.nwith = zeros(Int, nlds.nσ+nlds.nρ)
    nlds.nused = zeros(Int, nlds.nσ+nlds.nρ)
    nlds.mycut = [mycut_d; mycut_e]
    nlds.trust = nothing
  end
end

function load!{S}(nlds::NLDS{S})
  if !nlds.loaded
    bigA = nlds.W
    if nlds.nθ > 0
      bigA = [bigA spzeros(size(bigA, 1), nlds.nθ)]
    end
    computecuts!(nlds)
    bigA = [bigA; get(nlds.cuts_DE)]

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
    MathProgBase.optimize!(nlds.model)
    status = MathProgBase.status(nlds.model)
    objval = MathProgBase.getobjval(nlds.model)
    primal = MathProgBase.getsolution(nlds.model)
    dual   = mygetdual(nlds.model)
    x = primal[1:end-nlds.nθ]
    θ = primal[end-nlds.nθ+1:end]
    if !isempty(nlds.nwith)
      nlds.nwith += 1
      nlds.nused[dual[nlds.nπ+1:end] .> 1e-6] += 1
      nlds.trust = nothing # need to be recomputed
    end
    π = dual[nlds.πs]
    σ = dual[nlds.σs]
    ρ = dual[nlds.ρs]
    cuts_d = get(nlds.cuts_de)[nlds.σs-nlds.nπ]
    cuts_e = get(nlds.cuts_de)[nlds.ρs-nlds.nπ]
    nlds.prevsol = NLDSSolution(status, objval, x, θ, vec(π' * nlds.T), vecdot(π, nlds.h), vecdot(σ, cuts_d), vecdot(ρ, cuts_e))
  end
end

function getsolution(nlds::NLDS)
  solve!(nlds)
  get(nlds.prevsol)
end
