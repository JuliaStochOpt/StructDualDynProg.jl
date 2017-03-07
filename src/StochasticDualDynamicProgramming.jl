__precompile__()

module StochasticDualDynamicProgramming

using DocStringExtensions

using CutPruners

using JuMP
using StructJuMP

import Base.show, Base.isless

# Utils
include("mycomp.jl")
include("stats.jl")

# Abstract components
# Stopping Criterion
include("stopcrit.jl")
# NLDS Model
include("cutstore.jl")
include("solver.jl")
include("nlds.jl")
# SDDP Tree
include("node.jl")
include("sddptree.jl")

# SDDP algorithm
include("sddp.jl")

# Wait and See value
include("waitandsee.jl")

# StructJuMP interface
include("interface.jl")

end # module
