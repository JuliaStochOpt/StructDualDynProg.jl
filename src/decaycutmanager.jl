export DecayCutManager

type DecayCutManager{S} <: AbstractCutManager{S}
  # used to generate cuts
  cuts_DE::Nullable{AbstractMatrix{S}}
  cuts_de::Nullable{AbstractVector{S}}

  nσ::Int
  nρ::Int
  σs::Vector{Int}
  ρs::Vector{Int}

  maxncuts::Int

  trust::Vector{Float64}

  λ::Float64
  newcuttrust::Float64
  mycutbonus::Float64

  function DecayCutManager(maxncuts::Int, λ=0.9, newcuttrust=0.8, mycutbonus=1)#newcuttrust=(1/(1/0.9-1))/2, mycutbonus=(1/(1/0.9-1))/2)
    new(nothing, nothing, 0, 0, Int[], Int[], maxncuts, Float64[], λ, newcuttrust, mycutbonus)
  end
end

DecayCutManager(maxncuts::Int, λ=0.9, newcuttrust=(1/(1/0.9-1))/2, mycutbonus=(1/(1/0.9-1))/2) = DecayCutManager{Float64}(maxncuts, λ, newcuttrust, mycutbonus)

function clone{S}(man::DecayCutManager{S})
  DecayCutManager{S}(man.maxncuts, man.λ, man.newcuttrust, man.mycutbonus)
end

# COMPARISON
function updatestats!(man::DecayCutManager, σρ)
  if ncuts(man) > 0
    man.trust *= man.λ
    man.trust[σρ .> 1e-6] += 1
  end
end

function initialtrust(man::DecayCutManager, mycut)
  if mycut
    man.newcuttrust + man.mycutbonus
  else
    man.newcuttrust
  end
end

function isbetter(man::AbstractCutManager, i::Int, mycut::Bool)
  if mycut
    # If the cut has been generated, that means it is usefull
    false
  else
    # The new cut has initial trust initialtrust(man, false)
    # but it is a bit disadvantaged since it is new so
    # as we advantage the new cut if mycut == true,
    # we advantage this cut by taking initialtrust(man, true)
    # with true instead of false
    man.trust[i] > initialtrust(man, true)
  end
end
