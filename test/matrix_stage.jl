if !isdefined(@__MODULE__, :lp_solvers)
    include("solvers.jl")
end
import SparseArrays

@testset "Matrix stages with a transition map and $cutgen" for cutgen in
    (StructProg.MultiCutGenerator(), StructProg.AvgCutGenerator())
    pruner = AvgCutPruningAlgo(-1)
    root = StructProg.NLDS{Float64}(
        ones(1, 1), [2.0], zeros(1, 0), [(:Zero, [1])], [(:Free, [1])],
        [1.0], GLPK.Optimizer, pruner,
    )
    child = StructProg.NLDS{Float64}(
        SparseArrays.sparse(ones(1, 1)), [0.0], fill(-1.0, 1, 1), [(:NonPos, [1])],
        [(:NonNeg, [1])], [1.0], GLPK.Optimizer, pruner, :AddImmediately,
    )
    sp = StructProg.StochasticProgram{Float64}(2)
    r = SOI.add_scenario_node!(sp, StructProg.NodeData(root, 0))
    c = SOI.add_scenario_node!(sp, StructProg.NodeData(child, 1))
    SOI.set!(sp, StructProg.CutGenerator(), r, cutgen)
    SOI.add_scenario_transition!(sp, r, c, 1.0, fill(3.0, 1, 1))
    result = SOI.optimize!(sp, SDDP.Algorithm(K=1), SOI.CutLimit(0) | SOI.IterLimit(5))
    @test SOI.last_result(result).lowerbound ≈ 8
    @test SOI.last_result(result).upperbound ≈ 8
    @test SOI.optimize!(sp, WaitAndSee.Algorithm(GLPK.Optimizer, -1)) == (8.0, 0.0)
end
