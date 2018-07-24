using LightGraphs

## Stochastic Program

"""
    AbstractStochasticProgram <: LightGraphs.AbstractGraph{Int}

Stochastic program instance
"""
abstract type AbstractStochasticProgram <: LightGraphs.AbstractGraph{Int} end

"""
    stochasticprogram(args...)

Creates a stochastic program from the arguments
"""
function stochasticprogram end

## State

# The type of states is `Int` for compatibility with LightGraphs' nodes

"""
    add_scenario_state!(sp::AbstractStochasticProgram, ...)

Add a new state to the stochastic program `sp` and returns it.
"""
function add_scenario_state! end

"""
    addcut!(sp::AbstractStochasticProgram, state, pool::SOI.AbstractSolutionPool, stats, ztol)

Add cut `cut` to the state `state` using the solution pool `pool` and the threshold `ztol` to determine if the cut is redundant.
The statistics are recorded in `stats`.
"""
function addcut! end

"""
    applycuts!(sp::AbstractStochasticProgram, state)

Apply cuts additions to the state `state`.
"""
function applycuts! end

## Transition

"""
    AbstractTransition <: LightGraphs.AbstractEdge{Int}

Transition between two states of the stochastic program
"""
abstract type AbstractTransition <: LightGraphs.AbstractEdge{Int} end

"""
    add_scenario_transition!(sp::AbstractStochasticProgram, ...)

Add a new transition to the stochastic program `sp` and returns it.
"""
function add_scenario_transition! end
