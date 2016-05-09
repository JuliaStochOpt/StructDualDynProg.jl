type CutStore{S}
  A::AbstractMatrix{S}
  b::AbstractVector{S}

  followers::Vector

  function CutStore(nvars)
    cs = new(spzeros(S, 0, nvars), spzeros(S, 0), Vector{Tuple{NLDS{S},Tuple{Symbol,Int64}}}(0))
    cs
  end
end

function addcut{S}(store::CutStore{S}, a::Vector{S}, β::S)
  store.A = [store.A; sparse(a')]
  store.b = [store.b; sparsevec([β])]

  for follower in store.followers
    notifynewcut(follower[1], a, β, follower[2])
  end
end

function addfollower(store::CutStore, nlds)
  push!(store.followers, nlds)
end
