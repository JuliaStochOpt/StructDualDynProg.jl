# TODO This module should be replaced by the package https://github.com/JuliaStochOpt/StochOptInterface.jl
module StochOptInterface

using DocStringExtensions

# Stochastic Program
include("stochprog.jl")
include("attributes.jl")

# Stats
include("stats.jl")
# Stopping Criterion
include("stopcrit.jl")
# Solution
include("solution.jl")

# Algorithm
include("algorithm.jl")

end
