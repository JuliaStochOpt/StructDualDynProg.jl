module StochasticDualDynamicProgramming

using JuMP
using StructJuMP

include("cutstore.jl")
include("nlds.jl")
include("node.jl")
include("interface.jl")

end # module
