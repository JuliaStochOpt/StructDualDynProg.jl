type CutStore{S}
  A::AbstractMatrix{S}
  b::AbstractVector{S}
  authors::Vector

  followers::Vector

  function CutStore(nvars)
    new(spzeros(S, 0, nvars), spzeros(S, 0), Vector{NLDS{S}}(0), Vector{Tuple{NLDS{S},Tuple{Symbol,Int64}}}(0))
  end
end

function addcut{S}(store::CutStore{S}, a::Vector{S}, β::S, author)
  store.A = [store.A; sparse(a')]
  store.b = [store.b; sparsevec([β])]
  push!(store.authors, author)

  for follower in store.followers
    notifynewcut(follower[1], a, β, follower[2], author)
  end
end

function addfollower(store::CutStore, nlds)
  push!(store.followers, nlds)
end
