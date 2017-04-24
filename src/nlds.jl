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
    # Number of feasibility cuts
    nσ::Int
    # Location of each feasibility cut in the constraint matrix
    σs::Vector{Int}
    # Number of Optimality Cuts for each θ_i
    nρ::Vector{Int}
    # Location of each optimality cut in the constraint matrix
    ρs::Vector{Vector{Int}}

    model
    loaded
    solved
    prevsol::Nullable{NLDSSolution}

    newcut::Symbol
    pruningalgo::AbstractCutPruningAlgo
    FCpruner::AbstractCutPruner
    OCpruners::Vector

    function NLDS(W::AbstractMatrix{S}, h::AbstractVector{S}, T::AbstractMatrix{S}, K, C, c::AbstractVector{S}, solver, pruningalgo::AbstractCutPruningAlgo, newcut::Symbol=:AddImmediately)
        nx = size(W, 2)
        nθ = 0
        nπ = length(h)
        localFC = CutStore{S}(nx)
        localOC = CutStore{S}(nx)
        FCpruner = CutPruner{nx, S}(pruningalgo, :≥)
        OCpruners = typeof(FCpruner)[]
        if false
            model = MathProgBase.ConicModel(solver)
        else
            model = MathProgBase.LinearQuadraticModel(solver)
        end
        nlds = new(W, h, T, K, C, c, S[], nothing, nothing, CutStore{S}[], CutStore{S}[], localFC, localOC, nothing, nothing, :NoOptimalityCut, nx, nθ, nπ, 1:nπ, 0, Int[], Int[], Vector{Int}[], model, false, false, nothing, newcut, pruningalgo, FCpruner, OCpruners)
        addfollower(localFC, (nlds, (:Feasibility, 0)))
        addfollower(localOC, (nlds, (:Optimality, 1)))
        nlds
    end
end

function (::Type{NLDS{S}}){S}(W::AbstractMatrix, h::AbstractVector, T::AbstractMatrix, K, C, c::AbstractVector, solver, pruningalgo::AbstractCutPruningAlgo, newcut::Symbol=:AddImmediately)
    NLDS{S}(AbstractMatrix{S}(W), AbstractVector{S}(h), AbstractMatrix{S}(T), K, C, AbstractVector{S}(c), solver, pruningalgo, newcut)
end

function setchildren!{S}(nlds::NLDS{S}, childFC, childOC, proba, cutmode, childT)
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
    nlds.nρ = zeros(Int, nlds.nθ)
    nlds.ρs = [Int[] for i in 1:nlds.nθ]
    # θ ≧ max β - ⟨a, x⟩
    # so we give :Max and set lazy_minus to true because it is "- ⟨a, x⟩" and not "+ ⟨a, x⟩"
    nlds.OCpruners = [CutPruner{nlds.nx, S}(nlds.pruningalgo, :Max, true) for i in 1:nlds.nθ]
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

function appendchildren!{S}(nlds::NLDS{S}, childFC, childOC, proba, childT)
    @assert length(childFC) == length(childOC)
    @assert length(proba) == length(nlds.childOC) + length(childOC) == length(nlds.childFC) + length(childFC)
    oldnθ = nlds.nθ
    if nlds.cutmode == :MultiCut
        nlds.proba = proba
        nlds.nθ = length(proba)
    end
    Δθ = nlds.nθ - oldnθ
    append!(nlds.nρ, zeros(Int, Δθ))
    append!(nlds.ρs, [Int[] for i in 1:Δθ])
    append!(nlds.OCpruners, [CutPruner{nlds.nx, S}(nlds.pruningalgo) for i in 1:Δθ])
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

function updatemaxncuts!(nlds::NLDS, maxncuts, i)
    if i == 0
        nlds.FCpruner.maxncuts = maxncuts
    else
        nlds.OCpruners[i].maxncuts = maxncuts
    end
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
    i = attrs[2]
    if i > 0 && !isnull(nlds.childT)
        A = A * get(nlds.childT)[i]
    end
    # No .=== :(
    mine = [authors[i] === nlds for i in 1:length(authors)]
    if isfc
        pruner = nlds.FCpruner
    else
        pruner = nlds.OCpruners[i]
    end
    ncur = ncuts(pruner)
    addstatus = addcuts!(pruner, A, b, mine)
    npushed = sum(addstatus .> ncur)
    @assert ncur + npushed == ncuts(pruner)
    cur = nlds.nσ + sum(nlds.nρ)
    if isfc
        nlds.nσ += npushed
        is = nlds.σs
    else
        nlds.nρ[i] += npushed
        is = nlds.ρs[i]
    end
    append!(is, zeros(Int, npushed))
    for j in 1:nnewcuts
        if addstatus[j] > ncur
            cur += 1
            is[addstatus[j]] = cur
        end
        if addstatus[j] > 0 && nlds.loaded
            if 1 <= addstatus[j] <= ncur || nlds.newcut == :InvalidateSolver
                # Need to update cut
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
    @assert nlds.newcut == :InvalidateSolver || cur == nlds.nσ + sum(nlds.nρ)
    checkconsistency(nlds)
end

function checkconsistency(nlds)
    @assert length(nlds.πs) == nlds.nπ
    @assert length(nlds.σs) == nlds.nσ
    for i in 1:nlds.nθ
        @assert length(nlds.ρs[i]) == nlds.nρ[i]
    end
    ρs = reduce(append!, Int[], nlds.ρs)
    @assert sort([nlds.πs; nlds.nπ + nlds.σs; nlds.nπ + ρs]) == collect(1:(nlds.nπ + nlds.nσ + sum(nlds.nρ)))
end

function getrhs{S}(nlds::NLDS{S})
    checkconsistency(nlds)
    bs = [nlds.h - nlds.T * nlds.x_a]
    Ks = [nlds.K]
    Kcut = []
    cur = 0
    b = S[]
    if !isempty(nlds.FCpruner)
        append!(b, nlds.FCpruner.b)
        nlds.nσ = ncuts(nlds.FCpruner)
        nlds.σs = cur + (1:nlds.nσ)
        cur += ncuts(nlds.FCpruner)
        push!(Kcut, (:NonPos, nlds.σs))
    end
    for i in 1:nlds.nθ
        if !isempty(nlds.OCpruners[i])
            append!(b, nlds.OCpruners[i].b)
            nlds.nρ[i] = ncuts(nlds.OCpruners[i])
            nlds.ρs[i] = cur + (1:nlds.nρ[i])
            cur += ncuts(nlds.OCpruners[i])
            push!(Kcut, (:NonPos, nlds.ρs[i]))
        end
    end
    if !isempty(Kcut)
        push!(bs, b)
        push!(Ks, Kcut)
    end
    checkconsistency(nlds)
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

function computecuts!(nlds::NLDS)
    nlds.nπ = length(nlds.h)
    nlds.πs = collect(1:nlds.nπ)
    cuts_D, cuts_d, mycut_d = getfeasibilitycuts(nlds)
    cuts_D = padfcut(nlds, cuts_D)
    cuts_E, cuts_e, mycut_e = getoptimalitycuts(nlds)
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

function getcutsDE{S}(nlds::NLDS{S})
    checkconsistency(nlds)
    nc = nlds.nσ + sum(nlds.nρ)
    A = spzeros(S, nc, nlds.nx + nlds.nθ)
    if !isempty(nlds.FCpruner)
        A[nlds.σs, 1:nlds.nx] = nlds.FCpruner.A
    end
    for i in 1:nlds.nθ
        if !isempty(nlds.OCpruners[i])
            A[nlds.ρs[i], 1:nlds.nx] = nlds.OCpruners[i].A
            A[nlds.ρs[i], nlds.nx + i] = 1
        end
    end
    A
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

        # Needs to be done before getcutsDE
        bs, Ks = getrhs(nlds)

        computecuts!(nlds) # FIXME what is the use of this ?
        A = getcutsDE(nlds)
        if !isnull(nlds.xuray_a)
            A = [A spzeros(size(A, 1), 1)]
        end
        bigA = [bigA; A]

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

        myload!(nlds.model, bigc, bigA, bs, Ks, bigC)
        nlds.loaded = true
    end
end

function checkstatus(status::Symbol)
    if status == :Error
        error("The solver reported an error")
    elseif status == :UserLimit
        error("The solver reached iteration limit or timed out")
    end
end

function solve!{S}(nlds::NLDS{S})
    load!(nlds)
    if !nlds.solved
        MathProgBase.optimize!(nlds.model)
        status = MathProgBase.status(nlds.model)
        checkstatus(status)
        sol = NLDSSolution(status, _getobjval(nlds.model))

        if status == :Unbounded
            uray = _getunboundedray(nlds.model)
            sol.xuray = uray[1:nlds.nx]
            sol.objvalxuray = dot(nlds.c, sol.xuray)
            sol.θuray = uray[nlds.nx+(1:nlds.nθ)]

            # See https://github.com/JuliaOpt/Gurobi.jl/issues/80
            c = MathProgBase.getobj(nlds.model)

            MathProgBase.setobj!(nlds.model, zeros(c))
            MathProgBase.optimize!(nlds.model)
            newstatus = MathProgBase.status(nlds.model)
            @assert newstatus != :Unbounded # Now the objective is 0
            checkstatus(newstatus)
            if newstatus == :Infeasible
                # We discard unbounded ray, infeasibility is more important
                status = newstatus
                sol.xuray = nothing
                sol.objvalxuray = nothing
                sol.θuray = nothing
            else
                @assert newstatus == :Optimal
            end
            MathProgBase.setobj!(nlds.model, c)
        end

        if status != :Unbounded
            # if infeasible dual + λ iray is dual feasible for any λ >= 0
            # here I just take iray to build the feasibility cut
            dual = _getdual(nlds.model)
            if length(dual) == 0
                error("Dual vector returned by the solver is empty while the status is $status")
            end
            @assert length(dual) == nlds.nπ + nlds.nσ + sum(nlds.nρ)

            π = dual[nlds.πs]
            σρ = @view dual[nlds.nπ+1:end]

            σ = σρ[nlds.σs]
            addusage!(nlds.FCpruner, σ)
            sol.σd = vecdot(σ, nlds.FCpruner.b)

            sol.ρe = zero(S)
            for i in 1:nlds.nθ
                ρ = σρ[nlds.ρs[i]]
                addusage!(nlds.OCpruners[i], ρ)
                sol.ρe += vecdot(ρ, nlds.OCpruners[i].b)
            end

            sol.πT = vec(π' * nlds.T)
            sol.πh = vecdot(π, nlds.h)
        end

        # Needs to be done after since if status is unbounded I do a resolve
        if status != :Infeasible
            primal = _getsolution(nlds.model)
            sol.x = primal[1:nlds.nx]
            addposition!(nlds.FCpruner, sol.x)
            for i in 1:nlds.nθ
                addposition!(nlds.OCpruners[i], sol.x)
            end
            sol.objvalx = dot(nlds.c, sol.x)
            sol.θ = primal[nlds.nx+(1:nlds.nθ)]
        end

        nlds.prevsol = sol
    end
end

function getsolution(nlds::NLDS)
    solve!(nlds)
    get(nlds.prevsol)
end
