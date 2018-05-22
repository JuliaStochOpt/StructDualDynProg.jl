# Attributes

"""
    AbstractStochasticProgramAttribute

Abstract supertype for attribute objects that can be used to set or get attributes (properties) of the stochastic program.
"""
abstract type AbstractStochasticProgramAttribute end

"""
    AbstractStateAttribute

Abstract supertype for attribute objects that can be used to set or get attributes (properties) of states in the stochastic program.
"""
abstract type AbstractStateAttribute end

"""
    AbstractTransitionAttribute

Abstract supertype for attribute objects that can be used to set or get attributes (properties) of transitions in the stochastic program.
"""
abstract type AbstractTransitionAttribute end

"""
    get(sp::AbstractStochasticProgram, attr::AbstractStochasticProgramAttribute)

Return an attribute `attr` of the stochastic program `sp`.

    get(sp::AbstractStochasticProgram, attr::AbstractStateAttribute, state)

Return an attribute `attr` of the state `state` in stochastic program `sp`.

    get(sp::AbstractStochasticProgram, attr::AbstractTransitionAttribute, tr::AbstractTransition)

Return an attribute `attr` of the transition `tr` in stochastic program `sp`.

### Examples

```julia
get(model, MasterState())
get(model, Solution(), state)
get(model, Probability(), tr)
```
"""
function get end

"""
    set!(sp::AbstractStochasticProgram, attr::AbstractStochasticProgramAttribute, value)

Assign `value` to the attribute `attr` of the stochastic program `sp`.

    set!(sp::AbstractStochasticProgram, attr::AbstractStateAttribute, state, value)

Assign `value` to the attribute `attr` of the state `state` in stochastic program `sp`.

    set!(sp::AbstractStochasticProgram, attr::AbstractTransitionAttribute, tr::AbstractTransition, value)

Assign `value` to the attribute `attr` of the transition `tr` in stochastic program `sp`.

### Examples

```julia
set!(model, CutGenerator(), state, cutgen)
set!(model, SourceSolution(), tr, sol)
```
"""
function set! end

## Stochastic Program attributes

"""
    MasterState <: AbstractStochasticProgramAttribute

The master state.
"""
struct MasterState <: AbstractStochasticProgramAttribute end

"""
    struct NumberOfPaths <: AbstractStochasticProgramAttribute
        length::Int
    end

The number of paths of length `length` starting from the master state.
"""
struct NumberOfPaths <: AbstractStochasticProgramAttribute
    length::Int
end
get(sp::AbstractStochasticProgram, nop::NumberOfPaths) = get(sp, NumberOfPathsFrom(nop.length), get(sp, MasterState()))

## State attributes

# Graph-related
"""
    IsLeaf <: AbstractStateAttribute

Whether the state has no outgoing transitions.
"""
struct IsLeaf <: AbstractStateAttribute end

"""
    struct NumberOfPathsFrom <: AbstractStateAttribute
        length::Int
    end

The number of paths of length `length` starting from the state.
"""
struct NumberOfPathsFrom <: AbstractStochasticProgramAttribute
    length::Int
end

# Optimization-related
"""
    Solution <: AbstractStateAttribute

The solution at the state.
"""
struct Solution <: AbstractStateAttribute end

"""
    Dimension <: AbstractStateAttribute

The number of variables of the stochastic program at the state (not including the auxiliary variables used for the objective value of its outgoing transitions.
"""
struct Dimension <: AbstractStateAttribute end

"""
    CutGenerator <: AbstractStateAttribute

The cut generator of the state.
"""
struct CutGenerator <: AbstractStateAttribute end

"""
    StateObjectiveValueBound <: AbstractStateAttribute

The current bound to the objective of the state.

### Examples

If the program at state `state` is bounded and the objective value of its outgoing transtitions is bounded too (e.g. `TransitionObjectiveValueBound` has been set to a finite value), `MOI.get(sp, state)` returns a finite value summing.
"""
struct StateObjectiveValueBound <: AbstractStateAttribute end

## Transition attributes

"""
    Probability <: AbstractTransitionAttribute

The probability of the transition.
"""
struct Probability <: AbstractTransitionAttribute end

"""
    SourceSolution <: AbstractTransitionAttribute

The solution of the source of a transition to be used by its destination.
"""
struct SourceSolution <: AbstractTransitionAttribute end

"""
    TransitionObjectiveValueBound <: AbstractStateAttribute

The current bound to the objective of the state.
"""
struct TransitionObjectiveValueBound <: AbstractStateAttribute end
