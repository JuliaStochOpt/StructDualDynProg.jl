@testset "Optimize Stock with $solver" for solver in lp_solvers
    numScen = 2
    C = 1
    P = 2
    d = [2, 3]
    m1 = StructuredModel(num_scenarios=numScen)

    @variable m1 x >= 0
    @objective(m1, Min, C * x)

    for ξ in 1:numScen
        m2 = StructuredModel(parent=m1, prob=1/2, id=ξ)
        @variable(m2, s >= 0)
        @constraints m2 begin
            s <= d[ξ]
            s <= x
        end
        @objective(m2, Max, P * s)
    end

    num_stages = 2
    cutmode = :AveragedCut
    K = 2
    pereiracoef = 0.1

    root = model2lattice(m1, num_stages, solver, AvgCutPruningAlgo(-1), cutmode)
    sol = SDDP(root, num_stages, K = K, stopcrit = Pereira(0.1) | IterLimit(10), verbose = 0)

    # K = 10 is a multiple of 2 so with ProbaPathSampler(true), the sampling is deterministic
    # therefor we can test for sol.attrs[:niter]
    @test sol.attrs[:niter] == 4
    @test sol.status == :Optimal
    @test sol.objval == -2.0
end
