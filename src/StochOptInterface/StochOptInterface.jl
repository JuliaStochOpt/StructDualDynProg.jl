# TODO This module should be replaced by the package https://github.com/JuliaStochOpt/StochOptInterface.jl
module StochOptInterface

using DocStringExtensions

# Utils
include("comp.jl")

# Stochastic Program
include("stochprog.jl")
include("attributes.jl")

# Stats
include("stats.jl")
# Stopping Criterion
include("stopcrit.jl")
# Cut Generator
include("cutgen.jl")
# Solution
include("solution.jl")

end
