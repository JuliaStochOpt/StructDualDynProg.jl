function testniter_3stages(niter, K, maxncuts, cutgen, detectlb)
    if isa(cutgen, StructProg.MultiCutGenerator)
        if K == -1
            if maxncuts == -1
                if detectlb
                    # == 14 with Gurobi
                    # == 15 on Mac OS, == 17 on Linux
                    @test 14 <= niter <= 17
                else
                    # 14 on Gurobi, 15 on CPLEX and 19 otherwise
                    @test 14 <= niter <= 19
                end
            else
                if detectlb
                    # down to 14 with Gurobi
                    # down to 15 and up to 19 on Mac OS, 17 <= ... <= 18 on Linux
                    @test 14 <= niter <= 19
                else
                    # between 16, 17 on Gurobi
                    # between 15, 16 on CPLEX
                    # up to 22 on 32 bits, <= 22 on 64 bits
                    @test 15 <= niter <= 22
                end
            end
        else
            if detectlb
                # down to 11 with CPLEX
                # up to 16 on Windows, up to 15 on Mac OS, <= 14 on Linux
                @test 11 <= niter <= 16
            else
                # 12 on CPLEX
                # 13, 14, 21, 24 on Gurobi
                # up to 20 on Windows, up to 19 on Mac OS, <= 18 on Linux
                @test 12 <= niter <= 24
            end
        end
    else
        if maxncuts == -1
            if K == -1
                # 22 on Gurobi
                # == 19 on Mac OS, == 21 on Linux
                @test 19 <= niter <= 22
            else
                if detectlb
                    # == 14 on Mac OS, == 16 on Linux
                    @test 14 <= niter <= 16
                else
                    # == 14 on Mac OS, == 17 on Linux
                    @test 14 <= niter <= 17
                end
            end
        else
            if K == -1
                if detectlb
                    # 27 on Gurobi
                    # down to 15 Mac OS, >= 19 on Linux
                    @test 15 <= niter <= 27
                else
                    # 27 on Gurobi
                    @test 20 <= niter <= 27
                end
            else
                if detectlb
                    # down to 13 Mac OS, >= 15 on Linux
                    @test 13 <= niter <= 19
                else
                    # down to 14 Mac OS, >= 17 on Linux
                    @test 14 <= niter <= 19
                end
            end
        end
    end
end
