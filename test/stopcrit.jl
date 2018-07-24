using TimerOutputs

mutable struct InvalidStoppingCriterion <: SOI.AbstractStoppingCriterion
end

@testset "Stopping Criterion" begin
    K = 42
    z_LB = 1
    z_UB = 1
    σ = 1
    stats = StructDualDynProg.SOI.Stats()

    stats.upperbound = z_UB
    stats.lowerbound = z_LB
    stats.npaths = K
    stats.σ_UB = σ

    totalstats = StructDualDynProg.SOI.Stats()
    totalstats.niterations = 8

    to = TimerOutput()
    for i in 1:8
        @timeit to StructDualDynProg.SOI.OCUTS_KEY begin end
    end

    @test_throws ErrorException SOI.stop(InvalidStoppingCriterion(), to, stats, stats)
    @test SOI.stop(SOI.AndStoppingCriterion(SOI.IterLimit(8), SOI.CutLimit(8)), to, stats, totalstats)
    totalstats.niterations = 7
    to = TimerOutput()
    for i in 1:2
        @timeit to StructDualDynProg.SOI.OCUTS_KEY begin end
    end
    for i in 1:6
        @timeit to StructDualDynProg.SOI.FCUTS_KEY begin end
    end
    @test !SOI.stop(SOI.AndStoppingCriterion(SOI.IterLimit(8), SOI.CutLimit(8)), to, stats, totalstats)
    to = TimerOutput()
    for i in 1:2
        @timeit to StructDualDynProg.SOI.OCUTS_KEY begin end
    end
    for i in 1:5
        @timeit to StructDualDynProg.SOI.FCUTS_KEY begin end
    end
    @test !SOI.stop(SOI.AndStoppingCriterion(SOI.IterLimit(8), SOI.CutLimit(8)), to, stats, totalstats)
end
