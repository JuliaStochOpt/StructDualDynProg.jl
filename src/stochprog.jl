using LightGraphs

export AbstractStochasticProgram, stochasticprogram
export isleaf, probability, cutgenerator, numberofpaths
export setprobability!, add_scenario_state!, add_scenario_transition!

"""
    AbstractStochasticProgram

Stochastic program instance
"""
abstract type AbstractStochasticProgram <: AbstractGraph{Int} end

"""
    stochasticprogram(args...)

Creates a stochastic program from the arguments
"""
function stochasticprogram end

"""
    getmaster(sp::AbstractStochasticProgram)

Returns the master node of stochastic problem `sp`.
"""
function getmaster end

"""
    probability(sp::AbstractStochasticProgram, edge)

Returns the probability to take the edge `edge` in the stochastic problem `sp`.
"""
function probability end

"""
    setchildx!(sp::AbstractStochasticProgram, node, child, sol)

Sets the parent solution of `child` as `sol`, the solution obtained at `node`.
"""
function setchildx! end

"""
    solve!(sp::AbstractStochasticProgram, node)

Solves the program at node `node` in `sp` and returns the solution.
"""
function solve! end

"""
    getobjectivebound(sp::AbstractStochasticProgram, node)

Gets the current bound to the objective of `node`.
"""
function getobjectivebound end

"""
    setθbound!(sp::AbstractStochasticProgram, node, child, θlb)

Sets the bounds to the objective of the child `child` of `node` to `θlb`.
"""
function setθbound! end

"""
    statedim(sp::AbstractStochasticProgram, node)

Returns the dimension of the state at `node`.
"""
function statedim end

"""
    isleaf(sp::AbstractStochasticProgram, node)

Returns whether `node` has no outgoing edge in `sp`.
"""
isleaf(sp::AbstractStochasticProgram, node) = iszero(outdegree(sp, node))

"""
    numberofpaths(sp::AbstractStochasticProgram, node, len::Integer)

Returns number of paths of length `len` starting at node `node` in `sp`. It should return 1 of `len` is zero.
"""
numberofpaths(sp::AbstractStochasticProgram, len) = numberofpaths(sp, getmaster(sp), len)

"""
    cutgenerator(sp::AbstractStochasticProgram, node)

Returns the cut generator of node `node` in the stochastic problem `sp`.
"""
function cutgenerator end

"""
    edgeid(sp::AbstractStochasticProgram, edge)

Returns the id of `edge` for `src(edge)`. This id should be between 1 and `outdegree(src(edge))`. In case of multiple cuts, the value `θ[edgeid(sp, edge)]` at `src(edge)` corresponds to the lower bound to the objective value of `dst(edge)`.
"""
function edgeid end

# Modification (Optional)

"""
    setprobability!(sp::AbstractStochasticProgram, edge, proba)

Sets the probability to take the edge `edge` in the stochastic problem `sp`.
"""
function setprobability! end

"""
    setcutgenerator!(sp::AbstractStochasticProgram, node, cutgen::AbstractOptimalityCutGenerator)

Sets the cut generator of node `node` to `cutgen` in the stochastic problem `sp`.
"""
function setcutgenerator! end
