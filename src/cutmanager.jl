function CutManager
  # used to generate cuts
  cuts_DE::Nullable{AbstractMatrix{S}}
  cuts_de::Nullable{AbstractVector{S}}

  nσ::Int
  nρ::Int
  σs::Vector{Int}
  ρs::Vector{Int}

  maxncuts::Integer

  nwith::Vector{Int}
  nused::Vector{Int}
  mycut::Vector{Bool}
  trust::Nullable{Vector{Float64}}

  newcuttrust::Float64
  mycutbonus::Float64
  bettercut::Nullable{Function}
end
