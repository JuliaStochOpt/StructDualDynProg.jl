function myveccat{S}(b::AbstractVector{S}, β::S)
  [b; β]
end
function myveccat{S}(b::AbstractSparseVector{S}, β::S)
  if β == zero(S)
    # If only homogeneous cuts, b stays sparse
    [b; sparsevec([β])]
  else
    # At the first non-homogeneous cut, b stops being sparse
    [b; β]
  end
end

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
    # FIXME shouldn't do sparse
    store.A = [store.A; sparse(a')]
    store.b = myveccat(store.b, β)
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
