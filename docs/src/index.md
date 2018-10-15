# StructDualDynProg.jl Documentation

This packages aims at providing an implementation of SDDP that is both efficient and modular/flexible.
It features the following:

* Support for unfeasible problem by generating a feasibility cut.
* Support for unbounded problem by using an unbounded ray.
* Support for a variety of cut pruning algorithm through the [CutPruners](https://github.com/JuliaPolyhedra/CutPruners.jl) package.
* Support for any linear or conic solvers available through [MathProgBase](https://github.com/JuliaOpt/MathProgBase.jl); see [JuliaOpt's webpage](http://www.juliaopt.org/) for a list.
* Support modeling the problem using the [StructJuMP modeling interface](https://github.com/StructJuMP/StructJuMP.jl).
* Support specifying the problem using a low-level interface. This is used for example by the [EntropicCone](https://github.com/blegat/EntropicCone.jl) package.

The packages is built on top of [StochOptInterface (SOI)](https://github.com/JuliastochOpt/StochOptInterface.jl/).
It implements an representation of stochastic programming implementing SOI in the `StructProg` submodule
and provides the following function for transforming a [StructJuMP](https://github.com/StructJuMP/StructJuMP.jl)
model into an instance of this representation:
```@docs
StructDualDynProg.StochOptInterface.stochasticprogram(m::StructDualDynProg.StructProg.JuMP.Model, num_stages, solver, pruningalgo::StructDualDynProg.StructProg.CutPruners.AbstractCutPruningAlgo, cutgen::StructProg.AbstractOptimalityCutGenerator, detectlb::Bool, newcut::Symbol)
```

This packages also provides an implemention of the SDDP algorithm that can be run on any stochastic
program implementing the SOI interface:
```@docs
SDDP.Algorithm
```

## Contents
```@contents
Pages = ["quickstart.md", "tutorial.md"]
Depth = 2
```
