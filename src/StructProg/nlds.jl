export NLDS, updatemaxncuts!

mutable struct Solution <: SOI.AbstractSolution
    status::Symbol
    objval
    objvalx
    objvalxuray # unbounded ray
    x
    xuray # unbounded ray
    θ
    θuray # unbounded ray
    cut_slope
    cut_intercept
    function Solution(status::Symbol, objval)
        new(status, objval, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing)
    end
end

SOI.feasibility_cut(sol::Solution) = (sol.cut_slope, sol.cut_intercept)
SOI.optimality_cut(sol::Solution) = (sol.cut_slope, sol.cut_intercept)

SOI.getstatus(sol::Solution) = sol.status

SOI.getobjectivevalue(sol::Solution) = sol.objval

SOI.getnodeobjectivevalue(sol::Solution) = sol.objvalx

SOI.getnodevalue(sol::Solution) = sol.x

function SOI.getbellmanvalue(sol::Solution, i)
    sol.θ[i]
end

# Nested L-Shaped Decomposition Subproblem (NLDS)

# Primal
# min c x
#     + θ                         -> AvgCutGenerator
#     + sum θ_i                   -> MultiCutGenerator
#     + λ c_a x_a                 -> parent unbounded
#     h -   Tx_a - Wx ∈ K         -> parent bounded
#     h - λ Tx_a - Wx ∈ K         -> parent unbounded
#     Dx >= d                     (feasibility cuts)
#     Ex + θ   >= e               -> AvgCutGenerator
#     Ex + θ_i >= e               -> MultiCutGenerator
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
mutable struct NLDS{S}
    stage::StructJuMP.StructuredModel
    model::JuMP.Model
    variables::Vector{JuMP.VariableRef}
    stage_cost::JuMP.AffExpr
    parent_constraints::Vector{JuMP.ConstraintRef}
    ray_step::Union{Nothing, JuMP.VariableRef}
    θ::Vector{JuMP.VariableRef}
    x_a::Vector{S}

    childFC::Vector{CutStore{S}}
    childOC::Vector{CutStore{S}}
    localFC::CutStore{S}
    localOC::CutStore{S}
    proba::Vector{Float64}
    θbounds::Dict{Int, JuMP.ConstraintRef}
    θfree::BitSet
    θlb::Vector{Float64}
    childT::Union{Nothing, Vector{AbstractMatrix{S}}}
    cutgen::AbstractOptimalityCutGenerator

    nx::Int
    σs::Vector{JuMP.ConstraintRef}
    ρs::Vector{Vector{JuMP.ConstraintRef}}

    prevsol::Union{Nothing, Solution}

    newcut::Symbol
    pruningalgo::AbstractCutPruningAlgo
    FCpruner::AbstractCutPruner
    OCpruners::Vector

    function NLDS{S}(stage, solver, nparent, pruningalgo, newcut) where S
        if !(newcut in (:AddImmediately, :InvalidateSolver))
            throw(ArgumentError("Unknown cut update mode $newcut"))
        end
        model = JuMP.Model(solver)
        parent = @variable(model, [1:nparent])
        parent_constraints = @constraint(model, [i=1:nparent], parent[i] == 0)
        variables, cost = add_stage!(model, stage, parent)
        @objective(model, Min, cost)
        nx = length(variables)
        localFC = CutStore{S}(nx)
        localOC = CutStore{S}(nx)
        FCpruner = CutPruner{nx, S}(pruningalgo, :≥)
        OCpruners = typeof(FCpruner)[]
        nlds = new{S}(stage, model, variables, cost, parent_constraints,
            nothing, JuMP.VariableRef[], zeros(S, nparent), CutStore{S}[],
            CutStore{S}[], localFC, localOC, Float64[], Dict{Int, JuMP.ConstraintRef}(), BitSet(), Float64[],
            nothing, AvgCutGenerator(), nx, JuMP.ConstraintRef[],
            Vector{JuMP.ConstraintRef}[], nothing, newcut, pruningalgo, FCpruner, OCpruners)
        addfollower(localFC, (nlds, (:Feasibility, 0)))
        addfollower(localOC, (nlds, (:Optimality, 0)))
        return nlds
    end
end

function add_childT!(nlds, childT)
    if childT === nothing
        @assert nlds.childT === nothing
    else
        # If there isn't any child yet, node.childT is null
        if nlds.childT === nothing
            nlds.childT = [childT]
        else
            push!(nlds.childT, childT)
        end
    end
end

function add_scenario_transition!(nlds::NLDS{S}, childFC, childOC, proba, childT) where {S}
    push!(nlds.proba, proba)
    push!(nlds.θlb, 0)
    push!(nlds.θfree, length(nlds.θlb))
    if length(nlds.θ) < nθ(nlds.cutgen, nlds.proba)
        push!(nlds.θ, @variable(nlds.model))
        push!(nlds.ρs, JuMP.ConstraintRef[])
        # θ ≧ max β - ⟨a, x⟩
        # so we give :Max and set lazy_minus to true because it is "- ⟨a, x⟩" and not "+ ⟨a, x⟩"
        push!(nlds.OCpruners, CutPruner{nlds.nx, S}(nlds.pruningalgo, :Max, true))
    end
    @assert length(childFC.b) == 0
    addfollower(childFC, (nlds, (:Feasibility, length(nlds.childFC)+1)))
    push!(nlds.childFC, childFC)
    @assert length(childOC.b) == 0
    addfollower(childOC, (nlds, (:Optimality, length(nlds.childOC)+1)))
    push!(nlds.childOC, childOC)
    add_childT!(nlds, childT)
    update_bellman!(nlds)
end

function setprobability!(nlds::NLDS, i, proba)
    nlds.proba[i] = proba
    update_bellman!(nlds)
end

function update_bellman!(nlds)
    for (i, θ) in enumerate(nlds.θ)
        averaged = nlds.cutgen isa AvgCutGenerator
        weight = averaged ? 1.0 : nlds.proba[i]
        JuMP.set_objective_coefficient(nlds.model, θ, weight)
        bounded = averaged ? isempty(nlds.θfree) : !(i in nlds.θfree)
        if bounded
            bound = averaged ? dot(nlds.proba, nlds.θlb) : nlds.θlb[i]
            if haskey(nlds.θbounds, i)
                JuMP.set_normalized_rhs(nlds.θbounds[i], bound)
            else
                nlds.θbounds[i] = @constraint(nlds.model, θ >= bound)
            end
        elseif haskey(nlds.θbounds, i)
            JuMP.delete(nlds.model, pop!(nlds.θbounds, i))
        end
    end
    nlds.prevsol = nothing
end

function getobjectivebound(nlds::NLDS)
    if !isempty(nlds.θ) && !isempty(nlds.θfree)
        return -Inf
    end
    bound = JuMP.constant(nlds.stage_cost)
    for (coef, var) in JuMP.linear_terms(nlds.stage_cost)
        if iszero(coef)
            continue
        elseif JuMP.is_fixed(var)
            bound += coef * JuMP.fix_value(var)
        elseif coef > 0 && JuMP.has_lower_bound(var)
            bound += coef * JuMP.lower_bound(var)
        elseif coef < 0 && JuMP.has_upper_bound(var)
            bound += coef * JuMP.upper_bound(var)
        else
            return -Inf
        end
    end
    return isempty(nlds.θ) ? bound : bound + dot(nlds.proba, nlds.θlb)
end

function setθbound!(nlds::NLDS, i, bound)
    if isfinite(bound)
        delete!(nlds.θfree, i)
        nlds.θlb[i] = bound
    else
        push!(nlds.θfree, i)
        nlds.θlb[i] = 0
    end
    update_bellman!(nlds)
end

function updatemaxncuts!(nlds::NLDS, maxncuts, i)
    if i == 0
        nlds.FCpruner.maxncuts = maxncuts
    else
        nlds.OCpruners[i].maxncuts = maxncuts
    end
end

function notifynewcuts(nlds::NLDS{S}, A::AbstractMatrix{S}, b::AbstractVector{S}, attrs, authors::Vector{NLDS{S}}) where S
    @assert attrs[1] in [:Feasibility, :Optimality]
    isfc = attrs[1] == :Feasibility
    i = attrs[2]
    if i > 0 && nlds.childT !== nothing
        A = A * nlds.childT[i]
    end
    if !isfc && i == 0
        i = 1  # The local averaged cut uses the single Bellman variable.
    end
    mine = map(author -> author === nlds, authors)
    if isfc
        pruner = nlds.FCpruner
    else
        pruner = nlds.OCpruners[i]
    end
    refs = isfc ? nlds.σs : nlds.ρs[i]
    slots = addcuts!(pruner, A, b, mine)
    resize!(refs, ncuts(pruner))
    for (j, slot) in enumerate(slots)
        if slot <= 0
            continue
        end
        if nlds.newcut == :InvalidateSolver
            MOI.Utilities.reset_optimizer(nlds.model)
        end
        if isassigned(refs, slot)
            JuMP.delete(nlds.model, refs[slot])
        end
        lhs = dot(A[j, :], nlds.variables)
        if !isfc
            lhs += nlds.θ[i]
        end
        refs[slot] = @constraint(nlds.model, lhs >= b[j])
        nlds.prevsol = nothing
    end
end

function setparentx(nlds::NLDS, x_a, ray, ray_cost)
    nlds.x_a = x_a
    if nlds.ray_step !== nothing
        JuMP.delete(nlds.model, nlds.ray_step)
        nlds.ray_step = nothing
    end
    if ray !== nothing && (!iszero(ray) || ray_cost < 0)
        step = @variable(nlds.model, lower_bound = 0)
        for (i, con) in enumerate(nlds.parent_constraints)
            JuMP.set_normalized_coefficient(con, step, -ray[i])
        end
        JuMP.set_objective_coefficient(nlds.model, step, ray_cost)
        nlds.ray_step = step
    end
    for (i, con) in enumerate(nlds.parent_constraints)
        JuMP.set_normalized_rhs(con, x_a[i])
    end
    nlds.prevsol = nothing
end

function read_solution!(nlds::NLDS{S}, sol) where S
    model = nlds.model
    status = JuMP.termination_status(model)
    if status == MOI.INFEASIBLE && JuMP.dual_status(model) == MOI.INFEASIBILITY_CERTIFICATE
        sol.status = :Infeasible
        sol.objval = Inf
        sol.xuray = sol.θuray = sol.objvalxuray = nothing
    elseif status == MOI.OPTIMAL
        sol.x = JuMP.value.(nlds.variables)
        sol.θ = JuMP.value.(nlds.θ)
        sol.objvalx = JuMP.value(nlds.stage_cost)
        addposition!(nlds.FCpruner, sol.x)
        for pruner in nlds.OCpruners
            addposition!(pruner, sol.x)
        end
        if sol.status == :Optimal
            sol.objval = JuMP.objective_value(model)
        end
    else
        error("Cannot construct a stage solution with solver status $status")
    end
    if sol.status != :Unbounded
        sol.cut_slope = -S[JuMP.dual(c) for c in nlds.parent_constraints]
        sol.cut_intercept = JuMP.dual_objective_value(model) + dot(sol.cut_slope, nlds.x_a)
        addusage!(nlds.FCpruner, JuMP.dual.(nlds.σs))
        for (pruner, refs) in zip(nlds.OCpruners, nlds.ρs)
            addusage!(pruner, JuMP.dual.(refs))
        end
    end
    return
end

function getsolution(nlds::NLDS{S}) where S
    if nlds.prevsol !== nothing
        return nlds.prevsol
    end
    model = nlds.model
    objective = JuMP.objective_function(model)
    sol = Solution(:Optimal, Inf)
    # Before the first cut, an unbounded Bellman variable gives a known ray.
    free = findfirst(eachindex(nlds.θ)) do i
        return isempty(nlds.ρs[i]) && !haskey(nlds.θbounds, i) &&
               JuMP.coefficient(objective, nlds.θ[i]) > 0
    end
    if free === nothing
        JuMP.optimize!(model)
        has_ray = JuMP.primal_status(model) == MOI.INFEASIBILITY_CERTIFICATE
        if has_ray
            sol.xuray = JuMP.value.(nlds.variables)
            sol.θuray = JuMP.value.(nlds.θ)
            sol.objvalxuray = JuMP.value(nlds.stage_cost) - JuMP.constant(nlds.stage_cost)
        end
    else
        has_ray = true
        sol.xuray = zeros(S, nlds.nx)
        sol.θuray = zeros(S, length(nlds.θ))
        sol.θuray[free] = -1
        sol.objvalxuray = zero(S)
    end
    if has_ray
        sol.status = :Unbounded
        sol.objval = -Inf
        # Recover a feasible base point for the recession direction.
        JuMP.set_objective_function(model, zero(JuMP.AffExpr))
        try
            JuMP.optimize!(model)
            read_solution!(nlds, sol)
        finally
            JuMP.set_objective_function(model, objective)
        end
    else
        read_solution!(nlds, sol)
    end
    nlds.prevsol = sol
    return sol
end

# Matrix interface: h - W*x - T*parent belongs to K, and x belongs to C.
function NLDS{S}(W::AbstractMatrix, h, T::AbstractMatrix, K, C, c, solver,
                 pruningalgo, newcut::Symbol=:AddImmediately) where S
    parent = StructuredModel(num_scenarios=1)
    p = @variable(parent, [1:size(T, 2)])
    stage = StructuredModel(parent=parent, id=1)
    x = @variable(stage, [1:size(W, 2)])
    for (cone, indices) in C
        for i in indices
            if cone == :NonNeg
                JuMP.set_lower_bound(x[i], 0)
            elseif cone == :NonPos
                JuMP.set_upper_bound(x[i], 0)
            elseif cone == :Zero
                JuMP.fix(x[i], 0)
            elseif cone != :Free
                throw(ArgumentError("Unsupported variable cone $cone"))
            end
        end
    end
    lhs = W * x
    if !isempty(p)
        lhs += T * p
    end
    for (cone, indices) in K
        for i in indices
            if cone == :NonNeg
                @constraint(stage, lhs[i] <= h[i])
            elseif cone == :NonPos
                @constraint(stage, lhs[i] >= h[i])
            elseif cone == :Zero
                @constraint(stage, lhs[i] == h[i])
            elseif cone != :Free
                throw(ArgumentError("Unsupported constraint cone $cone"))
            end
        end
    end
    @objective(stage, Min, dot(c, x))
    return NLDS{S}(stage, solver, length(p), pruningalgo, newcut)
end
