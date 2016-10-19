export AvgCutManager

type AvgCutManager{S} <: AbstractCutManager{S}
  # used to generate cuts
  cuts_DE::Nullable{AbstractMatrix{S}}
  cuts_de::Nullable{AbstractVector{S}}

  nσ::Int
  nρ::Int
  σs::Vector{Int}
  ρs::Vector{Int}

  maxncuts::Int

  nwith::Vector{Int}
  nused::Vector{Int}
  mycut::Vector{Bool}
  trust::Nullable{Vector{Float64}}

  newcuttrust::Float64
  mycutbonus::Float64

  function AvgCutManager(maxncuts::Int, newcuttrust=3/4, mycutbonus=1/4)
    new(nothing, nothing, 0, 0, Int[], Int[], maxncuts, Int[], Int[], Bool[], nothing, newcuttrust, mycutbonus)
  end
end

AvgCutManager(maxncuts::Int, newcuttrust=3/4, mycutbonus=1/4) = AvgCutManager{Float64}(maxncuts, newcuttrust, mycutbonus)

function clone{S}(man::AvgCutManager{S})
  AvgCutManager{S}(man.maxncuts, man.newcuttrust, man.mycutbonus)
end

function init!(man::AvgCutManager, mycut_d, mycut_e)
  man.nwith = zeros(Int, man.nσ+man.nρ)
  man.nused = zeros(Int, man.nσ+man.nρ)
  man.mycut = [mycut_d; mycut_e]
  man.trust = nothing
end

# COMPARISON
function updatestats!(man::AvgCutManager, σρ)
  if ncuts(man) > 0
    man.nwith += 1
    man.nused[σρ .> 1e-6] += 1
    man.trust = nothing # need to be recomputed
  end
end

function gettrustof(man::AvgCutManager, nwith, nused, mycut)
  (nwith == 0 ? man.newcuttrust : nused / nwith) + (mycut ? man.mycutbonus : 0)
end
function initialtrust(man::AvgCutManager, mycut)
  gettrustof(man, 0, 0, mycut)
end
function gettrust(man::AvgCutManager)
  if isnull(man.trust)
    trust = man.nused ./ man.nwith
    trust[man.nwith .== 0] = man.newcuttrust
    trust[man.mycut] += man.mycutbonus
    man.trust = trust
  end
  get(man.trust)
end

# CHANGE

function keeponly!(man::AvgCutManager, K::Vector{Int})
  man.nwith = man.nwith[K]
  man.nused = man.nused[K]
  man.mycut = man.mycut[K]
  man.trust = gettrust(man)[K]
end

function replacecuts!(man::AvgCutManager, js::Vector{Int}, mycut::Vector{Bool})
  man.nwith[js] = 0
  man.nused[js] = 0
  man.mycut[js] = mycut
  gettrust(man)[js] = initialtrusts(man, mycut)
end

function pushcuts!(man::AvgCutManager, mycut::Vector{Bool})
  append!(man.nwith, zeros(length(mycut)))
  append!(man.nused, zeros(length(mycut)))
  append!(man.mycut, mycut)
  if !isnull(man.trust)
    append!(get(man.trust), initialtrusts(man, mycut))
  end
end
