# Copy a stage into either a subproblem or a wait-and-see path model.
function add_stage!(model, stage, parent_variables)
    variables = [
        JuMP.add_variable(model, stage.variables[i], stage.varnames[i])
        for i in 1:stage.nextvaridx
    ]
    function variable(ref)
        if ref.model === stage
            return variables[ref.idx]
        elseif ref.model === stage.parent
            return parent_variables[ref.idx]
        end
        error("A stage may only reference its own variables and its parent's variables")
    end
    expression(ref::StructJuMP.StructuredVariableRef) = variable(ref)
    function expression(aff::JuMP.GenericAffExpr)
        result = JuMP.AffExpr(JuMP.constant(aff))
        for (coef, ref) in JuMP.linear_terms(aff)
            JuMP.add_to_expression!(result, coef, variable(ref))
        end
        return result
    end
    expression(func::AbstractVector) = map(expression, func)
    for i in 1:stage.nextconidx
        con = stage.constraints[i]
        JuMP.add_constraint(model, JuMP.build_constraint(error, expression(con.func), con.set), stage.connames[i])
    end
    cost = convert(JuMP.AffExpr, expression(JuMP.objective_function(stage)))
    if JuMP.objective_sense(stage) == MOI.MAX_SENSE
        cost = -cost
    end
    return variables, cost
end
