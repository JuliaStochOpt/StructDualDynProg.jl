include("solvers.jl")

# solver independent tests
include("comp.jl")
include("stopcrit.jl")
include("paths.jl")

# solver dependent tests
include("optimize_stock.jl")
include("hydro_thermal_scheduling.jl")
include("prob5.2.jl")
include("matrix_stage.jl")
