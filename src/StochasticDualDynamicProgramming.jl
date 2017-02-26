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
# Cut Manager
#include("cutmanager.jl")
#include("avgcutmanager.jl")
#include("decaycutmanager.jl")
# Stopping Criterion
include("stopcrit.jl")

# SDDP algorithm
include("cutstore.jl")
include("solver.jl")
include("nlds.jl")
include("node.jl")
include("sddp.jl")

# Wait and See value
include("waitandsee.jl")

# StructJuMP interface
include("interface.jl")

end # module
