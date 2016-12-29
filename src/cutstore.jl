function myveccat{S}(b::AbstractVector{S}, β::S, force=false)
    push!(b, β)
    b
end
function myveccat{S}(b::AbstractSparseVector{S}, β::S, force=false)
    if force || β == zero(S)
        # If only homogeneous cuts, b stays sparse
        [b; sparsevec([β])]
    else
        # At the first non-homogeneous cut, b stops being sparse
        [b; β]
    end
end

# see https://github.com/JuliaLang/julia/issues/16661
function mymatcat{S}(A::AbstractMatrix{S}, a::AbstractVector{S})
    [A; a']
end
function mymatcat{S}(A::AbstractSparseMatrix{S}, a::AbstractSparseVector{S})
    [A; sparse(a')]
end

type CutStore{S}
    A::AbstractMatrix{S}
    b::AbstractVector{S}
    authors::Vector
    Anew::AbstractMatrix{S}
    bnew::AbstractVector{S}
    authorsnew::Vector

    followers::Vector
    needstored::Vector{Bool}

    storecuts::Symbol

    function CutStore(nvars)
        new(spzeros(S, 0, nvars), spzeros(S, 0), NLDS{S}[], spzeros(S, 0, nvars), spzeros(S, 0), NLDS{S}[], Vector{Tuple{NLDS{S},Tuple{Symbol,Int64}}}(0), Vector{Bool}(0), :IfNeededElseDelete)
    end
end

function checksparseness(a::Vector)
    if true || countnz(a) * 2 < length(a)
        sparse(a)
    else
        a
    end
end

function addcut{S}(store::CutStore{S}, a::Vector{S}, β::S, author)
    a = checksparseness(a)
    store.Anew = mymatcat(store.Anew, a)
    store.bnew = myveccat(store.bnew, β)
    push!(store.authorsnew, author)
end

function apply!{S}(store::CutStore{S})
    if !isempty(store.bnew)
        if store.storecuts == :Yes || (store.storecuts != :No && reduce(|, false, store.needstored))
            store.A = [store.A; store.Anew]
            store.b = [store.b; store.bnew]
            append!(store.authors, store.authorsnew)
        end

        for follower in store.followers
            notifynewcuts(follower[1], store.Anew, store.bnew, follower[2], store.authorsnew)
        end

        store.Anew = spzeros(S, 0, size(store.A, 2))
        store.bnew = spzeros(S, 0)
        store.authorsnew = NLDS{S}[]
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
    if store.storecuts == :IfNeededElseDelete && !reduce(|, false, store.needstored)
        store.A = spzeros(S, 0, size(store.A, 2))
        store.b = spzeros(S, 0)
        store.authors = NLDS{S}[]
    end
end
