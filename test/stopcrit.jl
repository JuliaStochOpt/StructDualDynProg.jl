type InvalidStoppingCriterion <: AbstractStoppingCriterion
end

K = 42
z_LB = 1
z_UB = 1
σ = 1
stats = StochasticDualDynamicProgramming.SDDPStats()

stats.upperbound = z_UB
stats.niterations = 8
stats.nocuts = 8
stats.lowerbound = z_LB
stats.npaths = K
stats.σ_UB = σ

@test_throws ErrorException stop(InvalidStoppingCriterion(), stats)
@test stop(AndStoppingCriterion(IterLimit(8), CutLimit(8)), stats)
stats.niterations = 7
stats.nocuts = 2
stats.nfcuts = 6
@test !stop(AndStoppingCriterion(IterLimit(8), CutLimit(8)), stats)
stats.nocuts = 2
stats.nfcuts = 5
@test !stop(AndStoppingCriterion(IterLimit(8), CutLimit(8)), stats)
