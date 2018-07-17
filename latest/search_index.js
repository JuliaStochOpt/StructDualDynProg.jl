var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#StructDualDynProg.SDDP-Tuple{StructDualDynProg.AbstractStochasticProgram,Any}",
    "page": "Home",
    "title": "StructDualDynProg.SDDP",
    "category": "method",
    "text": "SDDP(sp, num_stages; K, stopcrit, verbose, pathsampler, ztol, stopatinf, mergepaths, forwardcuts, backwardcuts)\n\n\nRuns the SDDP algorithms on the stochastic program given by sp. The algorithm will do iterations until stopcrit decides to stop or when the root node is infeasible. In each iterations, K paths will be explored up to num_stages stages. The paths will be selected according to pathsampler and equivalent paths might be merged if their difference is smaller than ztol and mergepaths is true. The parameter ztol is also used to check whether a new cut is useful. When a scenario is infeasible and stopatinf is true then no other scenario with the same ancestor is run. Note that since the order in which the different scenarios is run is not deterministic, this might introduce nondeterminism even if the sampling is deterministic. By default, the cuts are added backward. However, if forwardcuts is set to true and backwardcuts is set to false the cuts are added forward.\n\n\n\n"
},

{
    "location": "index.html#StructDualDynProg.stochasticprogram-Tuple{JuMP.Model,Any,Any,CutPruners.AbstractCutPruningAlgo,StructDualDynProg.AbstractOptimalityCutGenerator,Bool,Symbol}",
    "page": "Home",
    "title": "StructDualDynProg.stochasticprogram",
    "category": "method",
    "text": "stochasticprogram(m::JuMP.Model, num_stages, solver,\n                  pruningalgo::CutPruners.AbstractCutPruningAlgo,\n                  cutgen::AbstractOptimalityCutGenerator=MultiCutGenerator(),\n                  detectlb::Bool=true, newcut::Symbol=:InvalidateSolver)\n\nCreates a StochasticProgram from a StructJuMP model. The former can then be used by the SDDP algorithm. The master problem is assumed to have model m and the scenarios are considered up to num_stages stages. The pruningalgo is as defined in CutPruners. If cutgen is MultiCutGenerator, one variable θ_i is created for each scenario. Otherwise, if cutgen is AveragedCutGenerator, only one variable θ is created and it represents the expected value of the objective value of the scenarios. If cutgen is NoOptimalityCut then no θ is created, only use this option if the objective of all models is zero except for the master model.\n\n\n\n"
},

{
    "location": "index.html#StructDualDynProg.jl-Documentation-1",
    "page": "Home",
    "title": "StructDualDynProg.jl Documentation",
    "category": "section",
    "text": "This packages aims at providing an implementation of SDDP that is both efficient and modular/flexible. It features the following:Support for unfeasible problem by generating a feasibility cut.\nSupport for unbounded problem by using an unbounded ray.\nSupport for a variety of cut pruning algorithm through the CutPruners package.\nSupport for any linear or conic solvers available through MathProgBase; see JuliaOpt\'s webpage for a list.\nSupport modeling the problem using the StructJuMP modeling interface.\nSupport specifying the problem using a low-level interface. This is used for example by the EntropicCone package.The SDDP algorithm can be run from any node of the lattice of problems using the following function:SDDP(sp::AbstractStochasticProgram, num_stages; K::Int, stopcrit::AbstractStoppingCriterion, verbose, pathsampler::AbstractPathSampler, ztol, stopatinf, mergepaths, forwardcuts, backwardcuts)This lattice can be built from a StructJuMP model using the following function:stochasticprogram(m::JuMP.Model, num_stages, solver, pruningalgo::CutPruners.AbstractCutPruningAlgo, cutgen::AbstractOptimalityCutGenerator, detectlb::Bool, newcut::Symbol)"
},

{
    "location": "index.html#Contents-1",
    "page": "Home",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [\"quickstart.md\", \"tutorial.md\", \"stopcrit.md\"]\nDepth = 2"
},

{
    "location": "quickstart.html#",
    "page": "Quick Start",
    "title": "Quick Start",
    "category": "page",
    "text": ""
},

{
    "location": "quickstart.html#A-first-example-:-Production-Planning-1",
    "page": "Quick Start",
    "title": "A first example : Production Planning",
    "category": "section",
    "text": "In this quick start guide, we show how to run the FAST quick start example using this package. We guide you through each step of the modeling separately. An IJulia notebook of this example can be found in the examples folder.We start by setting the different constantsconst num_stages = 2\nconst numScen = 2\nconst C = 1\nconst P = 2\nconst d = [2, 3]We now model the master problem using StructJuMP.using StructJuMP\nm1 = StructuredModel(num_scenarios=numScen)\n@variable(m1, x >= 0)\n@objective(m1, Min, C * x)For each of the two scenarios we need to create a StructJuMP model specifying that m1 is the parent and that the scenario has probability 1/2.for ξ in 1:numScen\n    m2 = StructuredModel(parent=m1, prob=1/2, id=ξ)\n    @variable(m2, s >= 0)\n    @constraints m2 begin\n        s <= d[ξ]\n        s <= x\n    end\n    @objective(m2, Max, P * s)\nendThis structured model need to be transformed into an appropriate structure to run SDDP on it. This is achieved by stochasticprogram:using GLPKMathProgInterface\nconst solver = GLPKMathProgInterface.GLPKSolverLP()\nusing CutPruners\nconst pruner = AvgCutPruningAlgo(-1)\nusing StructDualDynProg\nsp = stochasticprogram(m1, num_stages, solver, pruner)In this example, we have chosen the GLPK solver but you can use any LP solver listed in the table of the JuliaOpt\'s webpage.You can now run the sddp algorithm on it using SDDP:sol = SDDP(sp, num_stages, K = 2, stopcrit = Pereira(0.1) | IterLimit(10))We are using 2 forward paths per iteration and we stop either after 10 iterations or once the pereira criterion is satisfied with alpha = 01.We can verify that the algorithm have found the right value by inspecting the solution:@show sol.objval # sol.objval = -2.0"
},

{
    "location": "tutorial.html#",
    "page": "Tutorial",
    "title": "Tutorial",
    "category": "page",
    "text": ""
},

{
    "location": "tutorial.html#Hydro-Thermal-Scheduling-1",
    "page": "Tutorial",
    "title": "Hydro Thermal Scheduling",
    "category": "section",
    "text": "In this tutorial, we show how to run the FAST tutorial example using this package. The big difference between this example and the example quickstart example is that in this example we will model serial independence. There will be 5 stages and 2 scenarios per stages except for the first stage which has only one scenario. Each pair of scenario will have the same parent.We start by setting the constants:const num_stages = 5\nconst numScen = 2\nconst C = 5\nconst V = 8\nconst d = 6\nconst r = [2, 10]We now create a matrix to store all the variables of all the models. This allows us to use the variables of other models from a given model. We also create an array of the first model of each stage to give play the role of parent for the models of the next stage.using StructJuMP\nx = Matrix{JuMP.Variable}(num_stages, numScen)\ny = Matrix{JuMP.Variable}(num_stages, numScen)\np = Matrix{JuMP.Variable}(num_stages, numScen)\nmodels = Vector{JuMP.Model}(num_stages)Now, we create all the models. Note that each model declares that its parent is the first model (i.e. the model ξ == 1) of the previous stage. Hence if it is not the first model, it also declares that it has the same children than the first model of its stage. This is how serial independence is modeled in StructJuMP.for s in 1:num_stages\n    for ξ in 1:(s == 1 ? 1 : numScen) # for the first stage there is only 1 scenario\n        if s == 1\n            model = StructuredModel(num_scenarios=numScen)\n        else\n            model = StructuredModel(parent=models[s-1], prob=1/2, same_children_as=(ξ == 1 ? nothing : models[s]), id=ξ, num_scenarios=(s == num_stages ? 0 : numScen))\n        end\n        x[s, ξ] = @variable(model, lowerbound=0, upperbound=V)\n        y[s, ξ] = @variable(model, lowerbound=0)\n        p[s, ξ] = @variable(model, lowerbound=0)\n        if s > 1\n            @constraint(model, x[s, ξ] <= x[s-1, 1] + r[ξ] - y[s, ξ])\n        else\n            @constraint(model, x[s, ξ] <= mean(r) - y[s, ξ])\n        end\n        @constraint(model, p[s, ξ] + y[s, ξ] >= d)\n        @objective(model, Min, C * p[s, ξ])\n        # models[s] contains the first model only\n        if ξ == 1\n            models[s] = model\n        end\n    end\nendWe now create the lattice, note that the master problem is models[1].using GLPKMathProgInterface\nconst solver = GLPKMathProgInterface.GLPKSolverLP()\nusing CutPruners\nconst pruner = AvgCutPruningAlgo(-1)\nusing StructDualDynProg\nsp = stochasticprogram(models[1], num_stages, solver, pruner)The SDDP algorithm can now be run on the lattice:sol = SDDP(sp, num_stages, K = 16, stopcrit = Pereira(2., 0.5) | IterLimit(10))"
},

{
    "location": "stopcrit.html#",
    "page": "Stopping Criterion",
    "title": "Stopping Criterion",
    "category": "page",
    "text": ""
},

{
    "location": "stopcrit.html#StructDualDynProg.stop-Tuple{StructDualDynProg.AbstractStoppingCriterion,StructDualDynProg.AbstractSDDPStats,StructDualDynProg.AbstractSDDPStats}",
    "page": "Stopping Criterion",
    "title": "StructDualDynProg.stop",
    "category": "method",
    "text": "stop(s, stats, totalstats)\n\n\nReturns whether the SDDP algorithm should stop. If totalstats.niterations is 0, no iteration has already been done, otherwise, the niterationsth iteration has just finished. This iteration used stats.npaths paths and generated stats.nfcuts (resp. stats.nocuts) new feasibility (resp. optimality) cuts. The lower bound is now totalstats.lowerbound and the upper bound has mean totalstats.upperbound and variance totalstats.σ_UB.\n\n\n\n"
},

{
    "location": "stopcrit.html#StructDualDynProg.OrStoppingCriterion",
    "page": "Stopping Criterion",
    "title": "StructDualDynProg.OrStoppingCriterion",
    "category": "type",
    "text": "type OrStoppingCriterion <: StructDualDynProg.AbstractStoppingCriterion\n\nStops if lhs or rhs want to stop.\n\n\n\n"
},

{
    "location": "stopcrit.html#StructDualDynProg.AndStoppingCriterion",
    "page": "Stopping Criterion",
    "title": "StructDualDynProg.AndStoppingCriterion",
    "category": "type",
    "text": "type AndStoppingCriterion <: StructDualDynProg.AbstractStoppingCriterion\n\nStops if lhs and rhs want to stop.\n\n\n\n"
},

{
    "location": "stopcrit.html#StructDualDynProg.IterLimit",
    "page": "Stopping Criterion",
    "title": "StructDualDynProg.IterLimit",
    "category": "type",
    "text": "type IterLimit <: StructDualDynProg.AbstractStoppingCriterion\n\nStops if iter ≧ limit.\n\n\n\n"
},

{
    "location": "stopcrit.html#StructDualDynProg.CutLimit",
    "page": "Stopping Criterion",
    "title": "StructDualDynProg.CutLimit",
    "category": "type",
    "text": "type CutLimit <: StructDualDynProg.AbstractStoppingCriterion\n\nStops if there was less than or equal to limit cuts added in the iteration. For instance, CutLimit(0) stops when there are no cuts added.\n\n\n\n"
},

{
    "location": "stopcrit.html#StructDualDynProg.Pereira",
    "page": "Stopping Criterion",
    "title": "StructDualDynProg.Pereira",
    "category": "type",
    "text": "type Pereira <: StructDualDynProg.AbstractStoppingCriterion\n\nStops if z_UB - α * σ/√K - tol < z_LB < z_UB + α * σ/√K + tol and σ / √K > β * max(1, |z_LB|))\n\n\n\n"
},

{
    "location": "stopcrit.html#Stopping-Criterion-1",
    "page": "Stopping Criterion",
    "title": "Stopping Criterion",
    "category": "section",
    "text": "stop(s::AbstractStoppingCriterion, stats::AbstractSDDPStats, totalstats::AbstractSDDPStats)\nOrStoppingCriterion\nAndStoppingCriterion\nIterLimit\nCutLimit\nPereira"
},

]}
