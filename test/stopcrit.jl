using TimerOutputs

mutable struct InvalidStoppingCriterion <: SOI.AbstractStoppingCriterion
end

@testset "Stopping Criterion" begin
    #K = 42
    z_LB = 1
    z_UB = 1
    σ = 1

    result = SOI.Result()
    result.upperbound = z_UB
    result.lowerbound = z_LB
    #result.npaths = K
    result.σ_UB = σ

    info = SOI.Info()
    @timeit info.timer "iteration 8" begin
        for i in 1:8
            push!(info.results, result)
            @timeit info.timer SOI.OCUTS_KEY begin end
        end
    end

    @test_throws ErrorException SOI.stop(InvalidStoppingCriterion(), info)
    @test SOI.stop(SOI.AndStoppingCriterion(SOI.IterLimit(8), SOI.CutLimit(8)), info)

    info = SOI.Info()
    for i in 1:7
        push!(info.results, result)
    end
    @timeit info.timer "iteration 7" begin
        for i in 1:2
            @timeit info.timer SOI.OCUTS_KEY begin end
        end
        for i in 1:5
            @timeit info.timer SOI.FCUTS_KEY begin end
        end
    end
    @test !SOI.stop(SOI.AndStoppingCriterion(SOI.IterLimit(8), SOI.CutLimit(8)), info)

    @timeit info.timer "iteration 7" begin
        @timeit info.timer SOI.FCUTS_KEY begin end
    end
    @test !SOI.stop(SOI.AndStoppingCriterion(SOI.IterLimit(8), SOI.CutLimit(8)), info)
end
