__precompile__()

module StochasticDualDynamicProgramming

using JuMP
using StructJuMP
import Base.show, Base.isless

include("mycomp.jl")

include("cutmanager.jl")
include("avgcutmanager.jl")
include("decaycutmanager.jl")

include("cutstore.jl")
include("solver.jl")
include("nlds.jl")
include("node.jl")
include("sddp.jl")
include("waitandsee.jl")
include("interface.jl")

end # module
