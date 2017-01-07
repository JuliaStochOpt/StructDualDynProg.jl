#using ECOS
#solver = ECOS.ECOSSolver(verbose=false)
#using Clp
#solver = Clp.ClpSolver()
#using Gurobi
#solver = Gurobi.GurobiSolver(OutputFlag=0)
using GLPKMathProgInterface
solver = GLPKSolverLP()
