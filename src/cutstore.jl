type CutStore{T}
  A::Matrix{T}
  b::Vector{T}

  followers::Vector #{NLDS{T}}

  function CutStore(nvars)
    cs = new(Matrix{T}(0, nvars), Vector{T}(0), Vector{NLDS{T}}(0))
    cs
  end
end

function addcut{T}(store::CutStore{T}, a::Vector{T}, β::T)
  store.A = [store.A; a']
  store.b = [store.b; β]

  for follower in store.followers
    notifynewcut(follower)
  end
end

function addfollower(store::CutStore, nlds)
  push!(store.followers, nlds)
end
