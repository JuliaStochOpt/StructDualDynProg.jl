@testset "2 stages with $solver" for solver in lp_solvers
    include("prob5.2_data.jl")

    numScen = 2
    M1 = StructuredModel(num_scenarios=numScen)

    @variable M1 x1[1:n] >= 0
    @variable M1 v1[1:n] >= 0
    @constraints M1 begin
        x1 .== v1
    end
    @objective(M1, Min, dot(ic, v1))

    for s in 1:numScen
        M2 = StructuredModel(parent=M1, prob=p2[s], id=s)
        @variable(M2, y2[1:n, 1:m] >= 0)
        @constraints M2 begin
            demand[j=1:m], sum(y2[:,j]) == D2[j,s]
            ylim[i=1:n], sum(y2[i,:]) <= x1[i]
        end
        @objective(M2, Min, dot(C, y2 * T))
    end

    # With detectlb = true and MultiCut
    # Iteration | Status    |     LB    |      UB    | #FC | #OC |
    #     1     | Optimal   |      0.0  |     Inf    |  2  |  0  |
    #     2     | Optimal   |  22338    | 1196514.84 |  0  |  2  |
    #     3     | Optimal   |  89893.68 |  676982.10 |  0  |  2  |
    #     4     | Optimal   | 182717.07 |  353268.97 |  0  |  2  |
    #     5     | Optimal   | 303678.32 |  366193.45 |  0  |  2  |
    #     6     | Optimal   | 311405.77 |  356197.41 |  0  |  2  |
    #     7     | Optimal   | 319203.65 |  344833.83 |  0  |  2  |
    #     8     | Optimal   | 334922.82 |  368285.45 |  0  |  2  |
    #     9     | Optimal   | 339047.44 |  349576.21 |  0  |  1  |
    #    10     | Optimal   | 340315.52 |  340315.52 |  0  |  0  |
    # With detectlb = true and AveragedCut
    #     1     | Optimal   |      0.0  |     Inf    | #FC | #OC |
    #     2     | Optimal   |  22338    | 1196514.84 |  2  |  0  |
    #     3     | Optimal   |  73386.05 |  660474.47 |  0  |  1  |
    #     4     | Optimal   | 171778.69 |  356387.31 |  0  |  1  |
    #     5     | Optimal   | 297813.60 |  391523.18 |  0  |  1  |
    #     6     | Optimal   | 309298.85 |  353098.55 |  0  |  1  |
    #     7     | Optimal   | 314734.72 |  344892.19 |  0  |  1  |
    #     8     | Optimal   | 331094.09 |  376701.57 |  0  |  1  |
    #     9     | Optimal   | 336691.25 |  346276.39 |  0  |  1  |
    #    10     | Optimal   | 337198.76 |  346727.58 |  0  |  1  |
    #    11     | Optimal   | 339235.89 |  344363.19 |  0  |  1  | -> Pereira stops here
    #    12     | Optimal   | 339779.51 |  340949.37 |  0  |  1  |
    #    13     | Optimal   | 340315.52 |  340315.52 |  0  |  0  |
    function testniter(niter, K, maxncuts, cutgen, detectlb)
        if isa(cutgen, StructProg.MultiCutGenerator)
            @test niter == 10
        else
            if maxncuts == -1
                if K == -1
                    if detectlb
                        @test niter == 13
                    else
                        # 10 on Gurobi/CPLEX, 14 otherwise
                        @test niter == 10 || niter == 14
                    end
                else
                    # 9 on Gurobi/CPLEX, 11 otherwise
                    @test niter == 9 || niter == 11 # Pereira stops earlier
                end
            else
                if K == -1
                    if detectlb
                        # up to 15 on Mac OS, <= 14 on Linux
                        @test 13 <= niter <= 15
                    else
                        # 10 on Gurobi/CPLEX, between 14 and 15 otherwise
                        @test 10 <= niter <= 15
                    end
                else
                    # 9 on Gurobi/CPLEX, between 11 and 12 otherwise
                    @test 9 <= niter <= 12 # Pereira stops earlier
                end
            end
        end
    end
    fulltest(M1, 2, 340315.52, [5085,1311,3919,854], 334687.754566, 15869.996575, testniter, solver)
    # root = model2lattice(M1, 2, solver, cutmode)
    # sol = SDDP(root, 2, cutmode, :All, verbose)
    #
    # v11value = sol.sol[1:4]
    # @test sol.status == :Optimal
    # @test abs(sol.objval - 340315.52) < 0.1
    # @test norm(v11value - [5085,1311,3919,854]) < 0.1
    # SDDPclear(M1)
end
