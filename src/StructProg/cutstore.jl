function _veccat(b::AbstractVector{S}, β::S, force=false) where S
    push!(b, β)
    b
end
function _veccat(b::AbstractSparseVector{S}, β::S, force=false) where S
    if force || β == zero(S)
        # If only homogeneous cuts, b stays sparse
        [b; sparsevec([β])]
    else
        # At the first non-homogeneous cut, b stops being sparse
        [b; β]
    end
end

# see https://github.com/JuliaLang/julia/issues/16661
function _matcat(A::AbstractMatrix{S}, a::AbstractVector{S}) where S
    [A; a']
end
function _matcat(A::AbstractSparseMatrix{S}, a::AbstractSparseVector{S}) where S
    [A; sparse(a')]
end

mutable struct CutStore{S}
    A::AbstractMatrix{S}
    b::AbstractVector{S}
    authors::Vector
    Anew::AbstractMatrix{S}
    bnew::AbstractVector{S}
    authorsnew::Vector

    followers::Vector
    needstored::Vector{Bool}

    storecuts::Symbol

    function CutStore{S}(nvars) where {S}
        # spzeros(S, 0) -> S[] : See julia#22225
        new{S}(spzeros(S, 0, nvars), S[], NLDS{S}[], spzeros(S, 0, nvars), S[], NLDS{S}[], Tuple{NLDS{S},Tuple{Symbol,Int}}[], Bool[], :IfNeededElseDelete)
    end
end

function checksparseness(a::AbstractVector)
    if true || countnz(a) * 2 < length(a)
        sparse(a)
    else
        a
    end
end

function addcut(store::CutStore{S}, a::AbstractVector{S}, β::S, author) where S
    a = checksparseness(a)
    store.Anew = _matcat(store.Anew, a)
    store.bnew = _veccat(store.bnew, β)
    push!(store.authorsnew, author)
end

function apply!(store::CutStore{S}) where S
    if !isempty(store.bnew)
        if store.storecuts == :Yes || (store.storecuts != :No && Compat.reduce(|, store.needstored, init=false))
            store.A = [store.A; store.Anew]
            store.b = [store.b; store.bnew]
            append!(store.authors, store.authorsnew)
        end

        for follower in store.followers
            notifynewcuts(follower[1], store.Anew, store.bnew, follower[2], store.authorsnew)
        end

        store.Anew = spzeros(S, 0, size(store.A, 2))
        store.bnew = S[] # See julia#22225
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
function noneedstored!(store::CutStore{S}, nlds) where S
    store.needstored[findfirst(f->f[1] === nlds, store.followers)] = false
    if store.storecuts == :IfNeededElseDelete && !Compat.reduce(|, store.needstored, init=false)
        store.A = spzeros(S, 0, size(store.A, 2))
        store.b = S[] # See julia#22225
        store.authors = NLDS{S}[]
    end
end
