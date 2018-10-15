## Hydro Thermal Scheduling

In this tutorial, we show how to run the [FAST tutorial example](https://web.stanford.edu/~lcambier/fast/tuto.php) using this package.
The big difference between this example and the quickstart example is that in this example we will model serial independence.
There will be 5 stages and 2 scenarios per stages except for the first stage which has only one scenario.
Each pair of scenario will have the same parent.
An IJulia notebook of this example can be found [in the examples folder](https://github.com/JuliaStochOpt/StructDualDynProg.jl/blob/master/examples/Hydro_Thermal_Scheduling.ipynb).

We start by setting the constants:
```julia
const num_stages = 5
const numScen = 2
const C = 5
const V = 8
const d = 6
const r = [2, 10]
```

We now create a matrix to store all the variables of all the models.
This allows us to use the variables of other models from a given model.
We also create an array of the first model of each stage to give play the role of parent for the models of the next stage.
```julia
using StructJuMP
x = Matrix{JuMP.Variable}(undef, num_stages, numScen)
y = Matrix{JuMP.Variable}(undef, num_stages, numScen)
p = Matrix{JuMP.Variable}(undef, num_stages, numScen)
models = Vector{JuMP.Model}(undef, num_stages)
```

Now, we create all the models.
Note that each model declares that its parent is the first model (i.e. the model `ξ == 1`) of the previous stage.
Hence if it is not the first model, it also declares that it has the same children than the first model of its stage.
This is how serial independence is modeled in [StructJuMP](https://github.com/StructJuMP/StructJuMP.jl).
```julia
using Statistics
for s in 1:num_stages
    for ξ in 1:(s == 1 ? 1 : numScen) # for the first stage there is only 1 scenario
        if s == 1
            model = StructuredModel(num_scenarios=numScen)
        else
            model = StructuredModel(parent=models[s-1], prob=1/2, same_children_as=(ξ == 1 ? nothing : models[s]), id=ξ, num_scenarios=(s == num_stages ? 0 : numScen))
        end
        x[s, ξ] = @variable(model, lowerbound=0, upperbound=V)
        y[s, ξ] = @variable(model, lowerbound=0)
        p[s, ξ] = @variable(model, lowerbound=0)
        if s > 1
            @constraint(model, x[s, ξ] <= x[s-1, 1] + r[ξ] - y[s, ξ])
        else
            @constraint(model, x[s, ξ] <= mean(r) - y[s, ξ])
        end
        @constraint(model, p[s, ξ] + y[s, ξ] >= d)
        @objective(model, Min, C * p[s, ξ])
        # models[s] contains the first model only
        if ξ == 1
            models[s] = model
        end
    end
end
```

We now create the lattice, note that the master problem is `models[1]`.
```julia
using GLPKMathProgInterface
const solver = GLPKMathProgInterface.GLPKSolverLP()
using CutPruners
const pruner = AvgCutPruningAlgo(-1)
using StructDualDynProg
const SOI = StructDualDynProg.SOI
sp = SOI.stochasticprogram(models[1], num_stages, solver, pruner)
```

The SDDP algorithm can now be run on the lattice:
```julia
algo = StructDualDynProg.SDDP.Algorithm(K = 16)
sol = SOI.optimize!(sp, algo, SOI.Pereira(2., 0.5) | SOI.IterLimit(10))
```
