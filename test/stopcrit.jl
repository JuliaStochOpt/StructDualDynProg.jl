mutable struct InvalidStoppingCriterion <: AbstractStoppingCriterion
end

@testset "Stopping Criterion" begin
    K = 42
    z_LB = 1
    z_UB = 1
    σ = 1
    stats = StructDualDynProg.SDDPStats()

    stats.upperbound = z_UB
    stats.nocuts = 8
    stats.lowerbound = z_LB
    stats.npaths = K
    stats.σ_UB = σ

    totalstats = StructDualDynProg.SDDPStats()
    totalstats.niterations = 8

    @test_throws ErrorException stop(InvalidStoppingCriterion(), stats, stats)
    @test stop(AndStoppingCriterion(IterLimit(8), CutLimit(8)), stats, totalstats)
    totalstats.niterations = 7
    stats.nocuts = 2
    stats.nfcuts = 6
    @test !stop(AndStoppingCriterion(IterLimit(8), CutLimit(8)), stats, totalstats)
    stats.nocuts = 2
    stats.nfcuts = 5
    @test !stop(AndStoppingCriterion(IterLimit(8), CutLimit(8)), stats, totalstats)
end
