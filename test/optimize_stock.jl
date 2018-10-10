# See https://web.stanford.edu/~lcambier/fast/demo.php
@testset "Optimize Stock with $solver" for solver in lp_solvers
    isclp(solver) && continue
    numScen = 2
    C = 1
    P = 2
    d = [2, 3]
    @testset "Optimize Stock directly" begin
        # With min ⟨[1, 1], [x, θ]⟩ s.t. no constraint with x >= 0, θ free
        # Clp is unable to return the unbounded ray which should be [0, -1]
        m1 = StructuredModel(num_scenarios=numScen)
        @variable(m1, x >= 0)
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
        cutmode = StructProg.AvgCutGenerator()
        K = 2
        pereiracoef = 0.1

        # detectlb should change anything here since m2 is min -2s and s is bounded above by x with is unknown so it is unbounded
        for forwardcuts in [false, true]
            for detectlb in [false, true]
                sp = SOI.stochasticprogram(m1, num_stages, solver, AvgCutPruningAlgo(-1), cutmode, detectlb)
                @test SOI.get(sp, SOI.NumberOfPaths(0)) == 1
                @test SOI.get(sp, SOI.NumberOfPaths(1)) == 2
                @test SOI.get(sp, SOI.NumberOfPaths(2)) == 2
                @test sprint(show, StructDualDynProg.StructProg.nodedata(sp, 1)) == "Node of 1 variables\n"
                algo = SDDP.Algorithm(K = K, forwardcuts = forwardcuts, backwardcuts = !forwardcuts)
                sol = SOI.optimize!(sp, algo, SOI.Pereira(0.1) | SOI.IterLimit(10), 0)
                sstats = sprint(show, SOI.last_result(sol))
                @test occursin("Lower Bound: ", sstats)
                @test occursin("Upper Bound: ", sstats)

                # K = 10 is a multiple of 2 so with ProbaPathSampler(true), the sampling is deterministic
                # therefore we can test for sol.attrs[:niter]
                #
                # Iteration | Status    |  LB  |  UB  | #OC |
                #     1     | Unbounded | -Inf |  0.0 |  1  |
                #     2     | Unbounded | -Inf | -5.0 |  1  |
                #     3     | Optimal   | -2.5 | -2.0 |  1  |
                #     4     | Optimal   | -2.0 | -2.0 |  0  |
                @test SOI.niterations(sol) == 4
                @test SOI.last_result(sol).status == :Optimal
                @test SOI.last_result(sol).lowerbound == -2.0
                StructProg.clear(m1)
            end
        end
    end

    # See https://github.com/blegat/StructDualDynProg.jl/issues/10
    if !isgrb(solver) && !iscpx(solver)
        @testset "Add a scenario latter" begin
            m1 = StructuredModel(num_scenarios=numScen)
            @variable(m1, x >= 0)
            @objective(m1, Min, C * x)

            m2 = StructuredModel(parent=m1, prob=1., id=1)
            @variable(m2, s >= 0)
            @constraints m2 begin
                s <= d[2]
                s <= x
            end
            @objective(m2, Max, P * s)

            num_stages = 2
            # Multicut wouldn't work since we are adding a node
            cutmode = StructProg.AvgCutGenerator()
            K = 2
            pereiracoef = 0.1

            sp = SOI.stochasticprogram(m1, num_stages, solver, AvgCutPruningAlgo(-1), cutmode)
            sol = SOI.optimize!(sp, SDDP.Algorithm(K=K), SOI.Pereira(0.1) | SOI.IterLimit(10), 0)
            @test SOI.niterations(sol) == 3
            @test SOI.last_result(sol).status == :Optimal
            @test SOI.last_result(sol).lowerbound == -3.0

            m3 = StructuredModel(parent=m1, id=2)
            @variable(m3, s >= 0)
            @constraints m3 begin
                s <= d[1]
                s <= x
            end
            @objective(m3, Max, P * s)

            root = 1
            newnode = StructProg.createnode(sp, m3, 1, 1, solver, root, AvgCutPruningAlgo(-1), cutmode)
            @test length(sp.out_transitions[1]) == 1
            SOI.set!(sp, SOI.Probability(), sp.out_transitions[1][1], 1/2)
            SOI.add_scenario_transition!(sp, root, newnode, 1/2)
            sol = SOI.optimize!(sp, SDDP.Algorithm(K=K), SOI.Pereira(0.1) | SOI.IterLimit(10), 0)
            # 2 on Mac OS and Windows, 3 otherwise
            @test 2 <= SOI.niterations(sol) <= 3
            @test SOI.last_result(sol).status == :Optimal
            @test SOI.last_result(sol).lowerbound == -2.0
        end
    end
end
