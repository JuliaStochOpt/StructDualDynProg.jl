using JuMP
using StructJuMP
using StochasticDualDynamicProgramming
using Base.Test

#using ECOS
#solver = ECOS.ECOSSolver(verbose=false)
using Clp
solver = Clp.ClpSolver()
#using Gurobi
#solver = Gurobi.GurobiSolver(OutputFlag=0)
#using GLPKMathProgInterface
#solver = GLPKSolverLP()

function fulltest(m, num_stages, objval, solval)
  for mccount in [:All, 42]
    for maxncuts in [-1, 9]
      for newcut in [:AddImmediately, :InvalidateSolver]
        for cutmode in [:MultiCut, :AveragedCut]
          root = model2lattice(m, num_stages, solver, cutmode, newcut, maxncuts)
          sol = SDDP(root, num_stages, mccount)

          v11value = sol.sol[1:4]
          @test sol.status == :Optimal
          @test abs(sol.objval - objval) / objval < (mccount == :All ? 1e-6 : .03)
          @test norm(v11value - solval) / norm(solval) < (mccount == :All ? 1e-6 : .3)
          SDDPclear(m)
        end
      end
    end
  end
end

include("prob5.2.jl")
include("prob5.2_3stages.jl")
include("prob5.2_3stages_serial.jl")
