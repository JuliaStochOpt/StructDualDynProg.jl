type CutStore{S}
  A::Matrix{S}
  b::Vector{S}

  followers::Vector

  function CutStore(nvars)
    cs = new(Matrix{S}(0, nvars), Vector{S}(0), Vector{Tuple{NLDS{S},Tuple{Symbol,Int64}}}(0))
    cs
  end
end

function addcut{S}(store::CutStore{S}, a::Vector{S}, β::S)
  store.A = [store.A; a']
  store.b = [store.b; β]

  for follower in store.followers
    notifynewcut(follower[1], a, β, follower[2])
  end
end

function addfollower(store::CutStore, nlds)
  push!(store.followers, nlds)
end
