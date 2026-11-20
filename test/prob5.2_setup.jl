if !isdefined(@__MODULE__, :lp_solvers)
    include("solvers.jl")
end

function fulltest(m, num_stages, objval, solval, ws, wsσ, testniter, solver)
    isclp(solver) && return
    for K in [-1, 40]
        for maxncuts in [-1, 7]
            for newcut in [:InvalidateSolver]#[:AddImmediately, :InvalidateSolver]
                for cutmode in [StructProg.MultiCutGenerator(), StructProg.AvgCutGenerator()]
                    for detectlb in [false, true]
                        for pruningalgo in [AvgCutPruningAlgo(maxncuts), DecayCutPruningAlgo(maxncuts), DeMatosPruningAlgo(maxncuts)]
                            isclp(solver) && K == -1 && maxncuts == 7 && cutmode == :MultiCut && !detectlb && isa(pruningalgo, DeMatosPruningAlgo) && continue
                            sp = SOI.stochasticprogram(m, num_stages, solver, pruningalgo, cutmode, detectlb, newcut)

                            μ, σ = SOI.optimize!(sp, WaitAndSee.Algorithm(solver, K))
                            @test abs(μ - ws) / ws < (K == -1 ? 1e-6 : .03)
                            @test abs(σ - wsσ) / wsσ <= (K == -1 ? 1e-6 : 1.)

                            if K == -1
                                stopcrit = SOI.CutLimit(0)
                            else
                                stopcrit = SOI.Pereira()
                            end
                            algo = SDDP.Algorithm(K = K, forwardcuts = true, backwardcuts = false)
                            sol = SOI.optimize!(sp, algo, stopcrit | SOI.IterLimit(42), 0)
                            if K == -1
                                testniter(SOI.niterations(sol), maxncuts, cutmode, detectlb)
                            else
                                # Sampled paths can reach the Pereira criterion at different iterations.
                                @test SOI.stop(stopcrit, sol)
                            end
                            @test SOI.last_result(sol).status == :Optimal
                            lb = SOI.last_result(sol).lowerbound
                            @test abs(lb - objval) / objval < (K == -1 ? 1e-6 : .03)
                            # FIXME the info is missing the SOI result
                            #v11value = sol.sol[1:4]
                            #@test norm(v11value - solval) / norm(solval) < (K == -1 ? 1e-6 : .3)

                            μ, σ = SOI.optimize!(sp, WaitAndSee.Algorithm(solver, K))
                            @test abs(μ - ws) / ws < (K == -1 ? 1e-6 : .03)
                            @test abs(σ - wsσ) / wsσ <= (K == -1 ? 1e-6 : 1.)

                        end
                    end
                end
            end
        end
    end
end
