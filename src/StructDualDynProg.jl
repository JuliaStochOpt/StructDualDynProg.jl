__precompile__()

module StructDualDynProg

using DocStringExtensions

using CutPruners

using JuMP
using StructJuMP

import Base.show, Base.isless

# Utils
include("comp.jl")

# Abstract components
# Stochastic Program
include("stochprog.jl")
# Stats
include("stats.jl")
# Stopping Criterion
include("stopcrit.jl")
# Cut Generator
include("cutgen.jl")
# Path Sampler
include("sampler.jl")
# Solution
include("solution.jl")

# SDDP algorithm on top of these abstract components
include("path.jl")
include("sddp.jl")

# Generic implementation of Stochastic Program
include("cutstore.jl")
include("solver.jl")
include("nlds.jl")
include("graph.jl")

# Wait and See value
include("waitandsee.jl")

# StructJuMP interface
include("interface.jl")

end # module
