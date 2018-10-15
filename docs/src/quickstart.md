## A first example : Production Planning

In this quick start guide, we show how to run the [FAST quick start example](https://web.stanford.edu/~lcambier/fast/demo.php) using this package.
We guide you through each step of the modeling separately.
An IJulia notebook of this example can be found [in the examples folder](https://github.com/JuliaStochOpt/StructDualDynProg.jl/blob/master/examples/Quick_Start.ipynb).

We start by setting the different constants
```julia
const num_stages = 2
const numScen = 2
const C = 1
const P = 2
const d = [2, 3]
```

We now model the master problem using [StructJuMP](https://github.com/StructJuMP/StructJuMP.jl).
```julia
using StructJuMP
m1 = StructuredModel(num_scenarios=numScen)
@variable(m1, x >= 0)
@objective(m1, Min, C * x)
```

For each of the two scenarios we need to create a [StructJuMP](https://github.com/StructJuMP/StructJuMP.jl) model specifying that `m1` is the parent and that the scenario has probability `1/2`.
```julia
for ξ in 1:numScen
    m2 = StructuredModel(parent=m1, prob=1/2, id=ξ)
    @variable(m2, s >= 0)
    @constraints m2 begin
        s <= d[ξ]
        s <= x
    end
    @objective(m2, Max, P * s)
end
```

This structured model need to be transformed into an appropriate structure to run SDDP on it.
This is achieved by [`StructDualDynProg.StochOptInterface.stochasticprogram`](@ref):
```julia
using GLPKMathProgInterface
const solver = GLPKMathProgInterface.GLPKSolverLP()
using CutPruners
const pruner = AvgCutPruningAlgo(-1)
using StructDualDynProg
using StochOptInterface
const SOI = StructDualDynProg.SOI
sp = SOI.stochasticprogram(m1, num_stages, solver, pruner)
```
In this example, we have chosen the [GLPK](https://github.com/JuliaOpt/GLPKMathProgInterface.jl/) solver but you can use any LP solver listed in the table of the [JuliaOpt's webpage](http://www.juliaopt.org/).

You can now run the [`SDDP.Algorithm`](@ref) on it:
```julia
algo = StructDualDynProg.SDDP.Algorithm(K = 2)
sol = SOI.optimize!(sp, algo, SOI.Pereira(0.1) | SOI.IterLimit(10))
```
We are using 2 forward paths per iteration and we stop either after 10 iterations or once the pereira criterion is satisfied with $\alpha = 0.1$.

We can verify that the algorithm have found the right value by inspecting the solution:
```julia
@show SOI.last_result(sol) # Lower and upper bounds are -2.0
```
