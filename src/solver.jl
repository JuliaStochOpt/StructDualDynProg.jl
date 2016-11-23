using MathProgBase

function mysetconstrB!(m::MathProgBase.AbstractConicModel, b, K)
    error("Not supported")
end

function getLPconstrbounds(bs, Ks)
    sumlen = sum(map(length, bs))
    lb = Vector{Float64}(sumlen)
    ub = Vector{Float64}(sumlen)
    offset = 0
    for i in 1:length(bs)
        b = bs[i]
        for (cone, idx) in Ks[i]
            #if !(cone in [:Zero, :NonPos, :NonNeg])
            #  error("This cone is not supported")
            #end
            offidx = offset+idx
            if cone == :Zero || cone == :NonPos
                lb[offidx] = b[idx]
            else
                lb[offidx] = -Inf
            end
            if cone == :Zero || cone == :NonNeg
                ub[offidx] = b[idx]
            else
                ub[offidx] = Inf
            end
        end
        offset += length(b)
    end
    lb, ub
end

function mysetconstrB!(m::MathProgBase.AbstractLinearQuadraticModel, bs, Ks)
    lb, ub = getLPconstrbounds(bs, Ks)
    MathProgBase.setconstrLB!(m, lb)
    MathProgBase.setconstrUB!(m, ub)
end

function myaddconstr!(m::MathProgBase.AbstractConicModel, idx, a, β, cone)
    error("Not supported")
end

function myaddconstr!(m::MathProgBase.AbstractLinearQuadraticModel, idx, a, β, cone)
    lb = -Inf
    ub = Inf
    if cone == :Zero || cone == :NonPos
        lb = β
    end
    if cone == :Zero || cone == :NonNeg
        ub = β
    end
    MathProgBase.addconstr!(m, idx, full(a), lb, ub)
end

function myload!(model::MathProgBase.AbstractConicModel, c, A, bs, Ks, C)
    MathProgBase.loadproblem!(nlds.model, c, A, reduce(vcat, Float64[], bs), reduce(vcat, [], Ks), C)
end

function myload!(model::MathProgBase.AbstractLinearQuadraticModel, c, A, bs, Ks, C)
    lb, ub = getLPconstrbounds(bs, Ks)
    l  = Vector{Float64}(size(A, 2))
    u  = Vector{Float64}(size(A, 2))
    for (cone, idx) in C
        if !(cone  in [:Free, :NonPos, :NonNeg])
            error("This cone is not supported")
        end
        if cone == :Free || cone == :NonPos
            l[idx] = -Inf
        else
            l[idx] = 0
        end
        if cone == :Free || cone == :NonNeg
            u[idx] = Inf
        else
            u[idx] = 0
        end
    end
    # TODO sparse objective triggers a bug in Gurobi solver, report it
    MathProgBase.loadproblem!(model, A, l, u, full(c), lb, ub, :Min)
end

function mygetdual(model::MathProgBase.AbstractConicModel)
    -MathProgBase.getdual(model)
end
function mygetdual(model::MathProgBase.AbstractLinearQuadraticModel)
    if MathProgBase.status(model) == :Infeasible
        MathProgBase.getinfeasibilityray(model)
    else
        MathProgBase.getconstrduals(model)
    end
end
