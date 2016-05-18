type CutStore{S}
  A::AbstractMatrix{S}
  b::AbstractVector{S}
  authors::Vector

  followers::Vector
  needstored::Vector{Bool}

  storecuts::Symbol

  function CutStore(nvars)
    new(spzeros(S, 0, nvars), spzeros(S, 0), Vector{NLDS{S}}(0), Vector{Tuple{NLDS{S},Tuple{Symbol,Int64}}}(0), Vector{Bool}(0), :IfNeededElseDelete)
  end
end

function addcut{S}(store::CutStore{S}, a::Vector{S}, β::S, author)
  if store.storecuts == :Yes || (store.storecuts != :No && reduce(|, false, store.needstored))
    store.A = [store.A; sparse(a')]
    store.b = [store.b; sparsevec([β])]
    push!(store.authors, author)
  end

  for follower in store.followers
    notifynewcut(follower[1], a, β, follower[2], author)
  end
end

function addfollower(store::CutStore, nlds)
  push!(store.followers, nlds)
  push!(store.needstored, true)
end

function needstored!(store::CutStore, nlds)
  store.needstored[findfirst(f->f[1] === nlds, store.followers)] = true
end
function noneedstored!{S}(store::CutStore{S}, nlds)
  store.needstored[findfirst(f->f[1] === nlds, store.followers)] = false
  if store.storecuts == :IfNeededElseDelete && reduce(|, false, store.needstored)
    store.A = spzeros(S, 0, size(store.A, 2))
    store.b = spzeros(S, 0)
  end
end
