export NLDS, updatemaxncuts!

# π, σ and ρ do not really make sense alone so only
# their product will T, h, d, e is stored
type NLDSSolution
    status::Symbol
    objval
    objvalx
    objvalxuray
    x
    xuray # unbouded ray
    θ
    θuray # unbouded ray
    πT
    πh
    σd
    ρe
    function NLDSSolution(status::Symbol, objval)
        new(status, objval, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing)
    end
end

function defaultbettercut(nwitha, nuseda, mycuta, nwithb, nusedb, mycutb)
    # If they are equal, then since we are asked whether "a" is striclty better than "b", we return false.
    # The stability of the sorting algorithm will do the rest.
    # For this reason, only strict inequalities are used ">", "<"

    # I use floating point arithmetic to avoid Int64 overflow
    if mycuta == mycutb
        if nwitha == 0 && nwithb == 0
            false
        elseif nwitha > 0 && nwithb == 0
            if mycutb
                false
            else
                nuseda / nwitha > 3/4
            end
        elseif nwitha == 0 && nwithb > 0
            if mycuta
                true
            else
                nusedb / nwithb < 3/4
            end
        else # nwitha > 0 && nwithb > 0
            nuseda / nwitha > nusedb / nwithb
        end
    elseif mycuta && !mycutb
        if nwitha == 0
            true
        elseif nwithb == 0 # nwitha > 0
            # if the ratio of "a" is larger than 1//2,
            # then we want to keep it over a new cut not created by this NLDS
            nuseda / nwitha > 1/2
        else # nwitha > 0 && nwithb > 0
            # "a" gets a 1/4 bonus trust because it was created by this NLDS
            nuseda / nwitha + 1/4 > nusedb / nwithb
        end
    else
        @assert mycutb && !mycuta
        !defaultbettercut(nwithb, nusedb, mycutb, nwitha, nuseda, mycuta)
    end
end

# Nested L-Shaped Decomposition Subproblem (NLDS)

# Primal
# min c x
#     + θ                         -> :AveragedCut
#     + sum θ_i                   -> :MultiCut
#     + λ c_a x_a                 -> parent unbounded
#     h -   Tx_a - Wx ∈ K         -> parent bounded
#     h - λ Tx_a - Wx ∈ K         -> parent unbounded
#     Dx >= d                     (feasibility cuts)
#     Ex + θ   >= e               -> :AveragedCut
#     Ex + θ_i >= e               -> :MultiCut
#     x in C
#     λ >= 0

# Note that we do not have θ >= 0

# Dual
# max π h + σ d + ρ e
#     - π Tx_a                    -> parent bounded
#     c - (π W + σ D + ρ E) ∈ C^* -> x dual
#     c_a x_a - π Tx_a >= 0       -> parent unbounded
#     π in K^*
#     1'ρ = 1 (or <= 1 if θ >= 0)
#     σ >= 0
#     ρ >= 0
type NLDS{S}
    W::AbstractMatrix{S}
    h::AbstractVector{S}
    T::AbstractMatrix{S}
    K
    C
    c::AbstractVector{S}

    # parent solution
    x_a::AbstractVector{S}
    xuray_a::Nullable{AbstractVector{S}}
    objvalxuray_a::Nullable{S}

    childFC::Vector{CutStore{S}}
    childOC::Vector{CutStore{S}}
    localFC::CutStore{S}
    localOC::CutStore{S}
    proba
    childT::Nullable{Vector{AbstractMatrix{S}}}
    cutmode::Symbol

    nx::Int
    nθ::Int
    nπ::Int
    πs::Vector{Int}

    model
    loaded
    solved
    prevsol::Nullable{NLDSSolution}

    newcut::Symbol
    pruner::AbstractCutPruner{S}

    function NLDS(W::AbstractMatrix{S}, h::AbstractVector{S}, T::AbstractMatrix{S}, K, C, c::AbstractVector{S}, solver, pruner::AbstractCutPruner{S}, newcut::Symbol=:AddImmediately)
        nx = size(W, 2)
        nθ = 0
        nπ = length(h)
        localFC = CutStore{S}(nx)
        localOC = CutStore{S}(nx)
        if false
            model = MathProgBase.ConicModel(solver)
        else
            model = MathProgBase.LinearQuadraticModel(solver)
        end
        nlds = new(W, h, T, K, C, c, S[], nothing, nothing, CutStore{S}[], CutStore{S}[], localFC, localOC, nothing, nothing, :NoOptimalityCut, nx, nθ, nπ, 1:nπ, model, false, false, nothing, newcut, CutPruner.clone(pruner))
        addfollower(localFC, (nlds, (:Feasibility, 0)))
        addfollower(localOC, (nlds, (:Optimality, 0)))
        nlds
    end
end

function (::Type{NLDS{S}}){S}(W::AbstractMatrix, h::AbstractVector, T::AbstractMatrix, K, C, c::AbstractVector, solver, pruner::AbstractCutPruner{S}, newcut::Symbol=:AddImmediately)
    NLDS{S}(AbstractMatrix{S}(W), AbstractVector{S}(h), AbstractMatrix{S}(T), K, C, AbstractVector{S}(c), solver, pruner, newcut)
end

function setchildren!(nlds::NLDS, childFC, childOC, proba, cutmode, childT)
    nlds.cutmode = cutmode
    @assert length(childFC) == length(childOC) == length(proba)
    if cutmode == :MultiCut
        nlds.proba = proba
        nlds.nθ = length(proba)
    elseif cutmode == :AveragedCut
        nlds.nθ = 1
    else
        nlds.nθ = 0
    end
    nlds.childFC = childFC
    for i in 1:length(childFC)
        addfollower(childFC[i], (nlds, (:Feasibility, i)))
    end
    nlds.childOC = childOC
    for i in 1:length(childOC)
        addfollower(childOC[i], (nlds, (:Optimality, i)))
    end
    nlds.childT = childT
end

function appendchildren!(nlds::NLDS, childFC, childOC, proba, childT)
    @assert length(childFC) == length(childOC)
    @assert length(proba) == length(nlds.childOC) + length(childOC) == length(nlds.childFC) + length(childFC)
    if nlds.cutmode == :MultiCut
        nlds.proba = proba
        nlds.nθ = length(proba)
    end
    for i in 1:length(childFC)
        @assert length(childFC[i].b) == 0
        addfollower(childFC[i], (nlds, (:Feasibility, length(nlds.childFC)+i)))
    end
    append!(nlds.childFC, childFC)
    for i in 1:length(childOC)
        @assert length(childFC[i].b) == 0
        addfollower(childOC[i], (nlds, (:Optimality, length(nlds.childOC)+i)))
    end
    append!(nlds.childOC, childOC)
    if childT === nothing
        @assert isnull(nlds.childT)
    else
        # If there isn't any child yet, nlds.childT is null
        if isnull(nlds.childT)
            nlds.childT = childT
        else
            append!(get(nlds.childT), childT)
        end
    end
end

function updatemaxncuts!(nlds::NLDS, maxncuts)
    nlds.pruner.maxncuts = maxncuts
end

# .=== doesn't work :(
function veceqeqeq(v::Vector, x)
    map(el->el === x, v)
end

function padfcut{S}(nlds::NLDS, D::AbstractMatrix{S})
    if nlds.nθ > 0
        [D spzeros(size(D, 1), nlds.nθ)]
    else
        D
    end
end

function getfeasibilitycuts(nlds::NLDS)
    function f(i)
        D = nlds.childFC[i].A
        if !isnull(nlds.childT)
            D = D * get(nlds.childT)[i]
        end
        D
    end
    cuts_D = reduce(vcat, nlds.localFC.A, map(i -> f(i), 1:length(nlds.childFC)))
    cuts_d = reduce(vcat, nlds.localFC.b, map(x -> x.b, nlds.childFC))
    mycut = reduce(vcat, veceqeqeq(nlds.localFC.authors, nlds), map(x -> veceqeqeq(x.authors, nlds), nlds.childFC))
    (cuts_D, cuts_d, mycut)
end

# see https://github.com/JuliaLang/julia/issues/16661
function padmultiocutaux{S}(nlds::NLDS, E::AbstractMatrix{S}, i, onesvec)
    nrows = size(E, 1)
    [E spzeros(S, nrows, i-1) onesvec spzeros(S, nrows, nlds.nθ-i)]
end
function padmultiocut{S}(nlds::NLDS, E::AbstractMatrix{S}, i)
    nrows = size(E, 1)
    padmultiocutaux(nlds, E, i, ones(S, nrows, 1))
end
function padmultiocut{S}(nlds::NLDS, E::AbstractSparseMatrix{S}, i)
    nrows = size(E, 1)
    padmultiocutaux(nlds, E, i, sparse(ones(S, nrows, 1)))
end

function padavgocut{S}(nlds::NLDS, E::AbstractMatrix{S})
    nrows = size(E, 1)
    [E ones(S, nrows)]
end
function padavgocut{S}(nlds::NLDS, E::AbstractSparseMatrix{S})
    nrows = size(E, 1)
    [E sparse(ones(S, nrows))]
end

function getoptimalitycuts{S}(nlds::NLDS{S})
    function f(i)
        E = nlds.childOC[i].A
        if !isnull(nlds.childT)
            E = E * get(nlds.childT)[i]
        end
        padmultiocut(nlds, E, i)
    end
    if nlds.nθ == 1
        cuts_E = padavgocut(nlds, nlds.localOC.A)
    else
        cuts_E = spzeros(S, 0, nlds.nx + nlds.nθ)
    end
    if nlds.nθ == length(nlds.childOC)
        cuts_E = reduce(vcat, cuts_E, map(f, 1:length(nlds.childOC)))
    end
    cuts_e = reduce(vcat, nlds.localOC.b, map(x -> x.b, nlds.childOC))
    mycut = reduce(vcat, veceqeqeq(nlds.localOC.authors, nlds), map(x -> veceqeqeq(x.authors, nlds), nlds.childOC))
    (cuts_E, cuts_e, mycut)
end


function notifynewcuts{S}(nlds::NLDS{S}, A::AbstractMatrix{S}, b::AbstractVector{S}, attrs, authors::Vector{NLDS{S}})
    @assert attrs[1] in [:Feasibility, :Optimality]
    isfc = attrs[1] == :Feasibility
    nnewcuts = size(A, 1)
    if isstarted(nlds.pruner)
        i = attrs[2]
        if i > 0 && !isnull(nlds.childT)
            A = A * get(nlds.childT)[i]
        end
        if isfc
            B = padfcut(nlds, A)
        else
            if i == 0
                B = padavgocut(nlds, A)
            else
                B = padmultiocut(nlds, A, i)
            end
        end
        # No .=== :(
        mine = [authors[i] === nlds for i in 1:length(authors)]
        addstatus = addcuts!(nlds.pruner, B, b, isfc, mine)
        for j in 1:nnewcuts
            if nlds.loaded
                if addstatus[j] != :Pushed || nlds.newcut == :InvalidateSolver
                    nlds.loaded = false
                    nlds.solved = false
                elseif nlds.newcut == :AddImmediately
                    idx = collect(1:nlds.nx)
                    a = A[j,:]
                    if !isfc
                        if i == 0
                            push!(idx, nlds.nx+1)
                        else
                            push!(idx, nlds.nx+i)
                        end
                        a = myveccat(a, one(S), true)
                    end
                    _addconstr!(nlds.model, idx, a, b[j], :NonPos)
                    nlds.solved = false
                else
                    error("Invalid newcut option $(nlds.newcut)")
                end
            end
        end
    end
end

function checkconsistency(nlds)
    @assert length(nlds.πs) == nlds.nπ
    @assert length(nlds.pruner.σs) == nlds.pruner.nσ
    @assert length(nlds.pruner.ρs) == nlds.pruner.nρ
    @assert sort([nlds.πs; nlds.nπ + nlds.pruner.σs; nlds.nπ + nlds.pruner.ρs]) == collect(1:(nlds.nπ + nlds.pruner.nσ + nlds.pruner.nρ))
end

function getrhs(nlds)
    bs = [nlds.h - nlds.T * nlds.x_a]
    Ks = [nlds.K]
    if ncuts(nlds.pruner) > 0
        push!(bs, get(nlds.pruner.cuts_de))
        Kcut = []
        if nlds.pruner.nσ > 0
            push!(Kcut, (:NonPos, nlds.pruner.σs))
        end
        if nlds.pruner.nρ > 0
            push!(Kcut, (:NonPos, nlds.pruner.ρs))
        end
        push!(Ks, Kcut)
    end
    bs, Ks
end

function setparentx(nlds::NLDS, x_a::AbstractVector, xuray_a, objvalxuray_a)
    nlds.x_a = x_a
    unbounded_a = xuray_a !== nothing
    if !isnull(nlds.xuray_a) || unbounded_a
        # FIXME do better when delvars!, ... are available in MPB
        nlds.loaded = false
        nlds.solved = false
    end
    if unbounded_a
        nlds.xuray_a = xuray_a
        nlds.objvalxuray_a = objvalxuray_a
    elseif !isnull(nlds.xuray_a)
        nlds.xuray_a = nothing
        nlds.objvalxuray_a = nothing
    end
    if nlds.loaded
        if unbounded_a
            nlds.loaded = false
        else
            bs, Ks = getrhs(nlds)
            mysetconstrB!(nlds.model, bs, Ks)
        end
        nlds.solved = false
    end
end

function computecuts!{S}(nlds::NLDS{S})
    if !isstarted(nlds.pruner)
        nlds.nπ = length(nlds.h)
        nlds.πs = collect(1:nlds.nπ)
        cuts_D, cuts_d, mycut_d = getfeasibilitycuts(nlds)
        cuts_D = padfcut(nlds, cuts_D)
        cuts_E, cuts_e, mycut_e = getoptimalitycuts(nlds)
        start!(nlds.pruner, cuts_D, cuts_E, cuts_d, cuts_e, mycut_d, mycut_e)
        # I will add the cuts in notifynewcuts! now
        noneedstored!(nlds.localFC, nlds)
        for store in nlds.childFC
            noneedstored!(store, nlds)
        end
        noneedstored!(nlds.localOC, nlds)
        for store in nlds.childOC
            noneedstored!(store, nlds)
        end
    end
end

function load!{S}(nlds::NLDS{S})
    if !nlds.loaded
        bigA = nlds.W
        if nlds.nθ > 0
            bigA = [bigA spzeros(size(bigA, 1), nlds.nθ)]
        end
        if !isnull(nlds.xuray_a)
            bigA = [bigA nlds.T * get(nlds.xuray_a)]
        end
        computecuts!(nlds)
        cuts_DE = get(nlds.pruner.cuts_DE)
        if !isnull(nlds.xuray_a)
            cuts_DE = [cuts_DE spzeros(size(cuts_DE, 1), 1)]
        end
        bigA = [bigA; cuts_DE]

        bigC = nlds.C
        bigc = nlds.c
        if nlds.nθ > 0
            bigC = [bigC; (:Free, collect(nlds.nx+(1:nlds.nθ)))]
            if nlds.proba === nothing
                @assert nlds.nθ == 1
                bigc = [bigc; 1]
            else
                bigc = [bigc; nlds.proba]
            end
        end
        if !isnull(nlds.xuray_a)
            bigC = [bigC; (:NonNeg, [nlds.nx+nlds.nθ+1])]
            bigc = [bigc; get(nlds.objvalxuray_a)]
        end

        bs, Ks = getrhs(nlds)

        myload!(nlds.model, bigc, bigA, bs, Ks, bigC)
        nlds.loaded = true
    end
end

function solve!(nlds::NLDS)
    load!(nlds)
    if !nlds.solved
        MathProgBase.optimize!(nlds.model)
        status = MathProgBase.status(nlds.model)
        if status == :Error
            error("The solver reported an error")
        elseif status == :UserLimit
            error("The solver reached iteration limit or timed out")
        else
            sol = NLDSSolution(status, _getobjval(nlds.model))
            if status != :Infeasible
                primal = _getsolution(nlds.model)
                sol.x = primal[1:nlds.nx]
                sol.objvalx = dot(nlds.c, sol.x)
                sol.θ = primal[nlds.nx+(1:nlds.nθ)]
            end
            if status == :Unbounded
                uray = _getunboundedray(nlds.model)
                sol.xuray = uray[1:nlds.nx]
                sol.objvalxuray = dot(nlds.c, sol.xuray)
                sol.θuray = uray[nlds.nx+(1:nlds.nθ)]
            else
                # if infeasible dual + λ iray is dual feasible for any λ >= 0
                # here I just take iray to build the feasibility cut
                dual = _getdual(nlds.model)
                if length(dual) == 0
                    error("Dual vector returned by the solver is empty while the status is $status")
                end
                @assert length(dual) == nlds.nπ + ncuts(nlds.pruner)

                π = dual[nlds.πs]
                σρ = dual[nlds.nπ+1:end]
                σ = σρ[nlds.pruner.σs]
                ρ = σρ[nlds.pruner.ρs]

                cuts_d = get(nlds.pruner.cuts_de)[nlds.pruner.σs]
                cuts_e = get(nlds.pruner.cuts_de)[nlds.pruner.ρs]

                sol.πT = vec(π' * nlds.T)
                sol.πh = vecdot(π, nlds.h)
                sol.σd = vecdot(σ, cuts_d)
                sol.ρe = vecdot(ρ, cuts_e)

                CutPruner.updatestats!(nlds.pruner, σρ)
            end

            nlds.prevsol = sol
        end
    end
end

function getsolution(nlds::NLDS)
    solve!(nlds)
    get(nlds.prevsol)
end
