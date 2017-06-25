function fulltest(m, num_stages, objval, solval, ws, wsσ, testniter, solver)
    isclp(solver) && return
    for K in [-1, 40]
        for maxncuts in [-1, 7]
            for newcut in [:InvalidateSolver]#[:AddImmediately, :InvalidateSolver]
                for cutmode in [MultiCutGenerator(), AvgCutGenerator()]
                    for detectlb in [false, true]
                        for pruningalgo in [AvgCutPruningAlgo(maxncuts), DecayCutPruningAlgo(maxncuts), DeMatosPruningAlgo(maxncuts)]
                            isclp(solver) && K == -1 && maxncuts == 7 && cutmode == :MultiCut && !detectlb && isa(pruningalgo, DeMatosPruningAlgo) && continue
                            root = model2lattice(m, num_stages, solver, pruningalgo, cutmode, detectlb, newcut)

                            μ, σ = waitandsee(root, num_stages, solver, K)
                            @test abs(μ - ws) / ws < (K == -1 ? 1e-6 : .03)
                            @test abs(σ - wsσ) / wsσ <= (K == -1 ? 1e-6 : 1.)

                            if K == -1
                                stopcrit = CutLimit(0)
                            else
                                stopcrit = Pereira()
                            end
                            stopcrit |= IterLimit(42)
                            sol = SDDP(root, num_stages, K = K, stopcrit = stopcrit, verbose = 0, forwardcuts = true, backwardcuts = false)
                            testniter(sol.attrs[:niter], K, maxncuts, cutmode, detectlb)
                            v11value = sol.sol[1:4]
                            @test sol.status == :Optimal
                            @test abs(sol.objval - objval) / objval < (K == -1 ? 1e-6 : .03)
                            @test norm(v11value - solval) / norm(solval) < (K == -1 ? 1e-6 : .3)

                            μ, σ = waitandsee(root, num_stages, solver, K)
                            @test abs(μ - ws) / ws < (K == -1 ? 1e-6 : .03)
                            @test abs(σ - wsσ) / wsσ <= (K == -1 ? 1e-6 : 1.)

                            SDDPclear(m)
                        end
                    end
                end
            end
        end
    end
end

@testset "Problem 5.2" begin
    include("prob5.2_2stages.jl")
    include("prob5.2_3stages.jl")
    include("prob5.2_3stages_serial.jl")
end
