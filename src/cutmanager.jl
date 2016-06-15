abstract AbstractCutManager{S}

function ncuts(man::AbstractCutManager)
  return length(get(man.cuts_de))
end

function isfeasibilitycut(man::AbstractCutManager, cut)
  if length(man.σs) < length(man.ρs)
    cut in man.σs
  else
    !(cut in man.ρs)
  end
end

function start!(man::AbstractCutManager, cuts_D, cuts_E, cuts_d, cuts_e, mycut_d, mycut_e)
  man.nσ = length(cuts_d)
  man.nρ = length(cuts_e)
  man.cuts_DE = [cuts_D; cuts_E]
  man.cuts_de = [cuts_d; cuts_e]
  man.σs = collect(1:man.nσ)
  man.ρs = collect(man.nσ+(1:man.nρ))
  init!(man, mycut_d, mycut_e)
end

function isstarted(man::AbstractCutManager)
  @assert isnull(man.cuts_DE) == isnull(man.cuts_de)
  !isnull(man.cuts_DE)
end

function choosecutstoremove(man::AbstractCutManager, num)
  # MergeSort is stable so in case of equality, the oldest cut loose
  # However PartialQuickSort is a lot faster

  trust = gettrust(man)
  if num == 1
    [indmin(trust)]                   # indmin selects the oldest cut in case of tie -> good :)
  else
    sortperm(trust, alg=PartialQuickSort(num))[1:num] # PartialQuickSort is unstable ->  bad :(
  end
end

isbetter(man::AbstractCutManager, i::Int, mycut::Bool) = gettrust(man)[i] > initialtrust(man, mycut)

# Add cut ax >= β
# If fc then it is a feasibility cut, otherwise it is an optimality cut
# If mycut then the cut has been added because of one of my trials
function addcut!{S}(man::AbstractCutManager{S}, a::AbstractVector{S}, β::S, isfc::Bool, mycut::Bool)

  if man.maxncuts != -1 && ncuts(man) >= man.maxncuts
    # Need to remove some cuts
    J = choosecutstoremove(man, ncuts(man) - man.maxncuts + 1)
    #if man.trustman.bettercut(0, 0, mycut, man.nwith[J[end]], man.nused[J[end]], man.mycut[J[end]])
    if !isbetter(man, J[end], mycut)
      j = J[end]
      get(man.cuts_DE)[j,:] = a
      get(man.cuts_de)[j] = β
      replacecut!(man, j, mycut)
      cutadded = true
      needupdate_σsρs = isfc $ isfeasibilitycut(man, j)
    else
      cutadded = false
      needupdate_σsρs = false
    end
    J = J[1:end-1]

    if length(J) > 1 || needupdate_σsρs
      keep = ones(Bool, ncuts(man))
      keep[J] = false
      K = find(keep)
      isσcut = zeros(Bool, ncuts(man))
      isσcut[man.σs] = true
      if cutadded
        isσcut[j] = isfc
      end
      isσcut = isσcut[K]
      man.σs = (1:length(isσcut))[isσcut]
      man.ρs = (1:length(isσcut))[!isσcut]
      man.nσ = length(man.σs)
      man.nρ = length(man.ρs)
    end

    if length(J) > 1
      man.cuts_DE = get(man.cuts_DE)[K,:]
      man.cuts_de = get(man.cuts_de)[K]
      keeponly!(man, K)
    end

    cutadded ? :Replaced : :Ignored
  else
    # Just append cut
    if isfc
      man.nσ += 1
      push!(man.σs, man.nσ + man.nρ)
    else
      man.nρ += 1
      push!(man.ρs, man.nσ + man.nρ)
    end
    man.cuts_DE = mymatcat(get(man.cuts_DE), a)
    man.cuts_de = myveccat(get(man.cuts_de), β)
    pushcut!(man, mycut)
    :Pushed
  end
end
