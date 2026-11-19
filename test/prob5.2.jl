using Test

@testset "Problem 5.2" begin
    include("prob5.2_2stages.jl")
    include("prob5.2_3stages.jl")
    include("prob5.2_3stages_serial.jl")
end
