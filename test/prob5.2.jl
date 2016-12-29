function fulltest(m, num_stages, objval, solval, ws, wsσ)
    for mccount in [-1, 42]
        for maxncuts in [-1, 7]
            for newcut in [:AddImmediately, :InvalidateSolver]
                for cutmode in [:MultiCut, :AveragedCut]
                    for cutmanager in [AvgCutManager(maxncuts), DecayCutManager(maxncuts)]
                        root = model2lattice(m, num_stages, solver, cutmanager, cutmode, newcut)

                        μ, σ = waitandsee(root, num_stages, solver, mccount)
                        @test abs(μ - ws) / ws < (mccount == -1 ? 1e-6 : .03)
                        @test abs(σ - wsσ) / wsσ <= (mccount == -1 ? 1e-6 : 1.)

                        sol = SDDP(root, num_stages, mccount=mccount, verbose=0, stopcrit=(x,y)->x>100)
                        v11value = sol.sol[1:4]
                        @test sol.status == :Optimal
                        @test abs(sol.objval - objval) / objval < (mccount == -1 ? 1e-6 : .03)
                        @test norm(v11value - solval) / norm(solval) < (mccount == -1 ? 1e-6 : .3)

                        μ, σ = waitandsee(root, num_stages, solver, mccount)
                        @test abs(μ - ws) / ws < (mccount == -1 ? 1e-6 : .03)
                        @test abs(σ - wsσ) / wsσ <= (mccount == -1 ? 1e-6 : 1.)

                        SDDPclear(m)
                    end
                end
            end
        end
    end
end

@testset "Problem 5.2" begin
    include("prob5.2_2stages.jl")
    #include("prob5.2_3stages.jl")
    #include("prob5.2_3stages_serial.jl")
end
