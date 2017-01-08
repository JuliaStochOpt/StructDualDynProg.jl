type InvalidStoppingCriterion <: AbstractStoppingCriterion
end

K = 42
z_LB = 1
z_UB = 1
σ = 1
@test_throws ErrorException stop(InvalidStoppingCriterion(), 1, 2, 3, K, z_LB, z_UB, σ)
@test stop(AndStoppingCriterion(IterLimit(8), CutLimit(8)), 8, 2, 6, K, z_LB, z_UB, σ)
@test !stop(AndStoppingCriterion(IterLimit(8), CutLimit(8)), 7, 2, 6, K, z_LB, z_UB, σ)
@test !stop(AndStoppingCriterion(IterLimit(8), CutLimit(8)), 7, 2, 5, K, z_LB, z_UB, σ)
