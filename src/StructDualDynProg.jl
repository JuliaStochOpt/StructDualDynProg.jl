__precompile__()

module StructDualDynProg

using Compat

using DocStringExtensions

using CutPruners

using JuMP
using StructJuMP

import Base.show, Base.isless

# Utils
include("comp.jl")
include("stats.jl")

# Abstract components
# Stopping Criterion
include("stopcrit.jl")
# Cut Generator
include("cutgen.jl")
# NLDS Model
include("cutstore.jl")
include("solver.jl")
include("nlds.jl")
# SDDP Graph
include("node.jl")
include("graph.jl")
# Path Sampler
include("sampler.jl")

# SDDP algorithm
include("path.jl")
include("sddp.jl")

# Wait and See value
include("waitandsee.jl")

# StructJuMP interface
include("interface.jl")

end # module
