# StochasticDualDynamicProgramming.jl Documentation

This packages aims at providing an implementation of SDDP that is both efficient and modular/flexible.
It features the following:

* Support for unfeasible problem by generating a feasibility cut.
* Support for unbounded problem by using an unbounded ray.
* Support for a variety of cut pruning algorithm through the [CutPruners](https://github.com/JuliaPolyhedra/CutPruners.jl) package.
* Support for any linear or conic solvers available through [MathProgBase](https://github.com/JuliaOpt/MathProgBase.jl); see [JuliaOpt's webpage](http://www.juliaopt.org/) for a list.
* Support modeling the problem using the [StructJuMP modeling interface](github.com/StructJuMP/StructJuMP.jl).
* Support specifying the problem using a low-level interface. This is used for example by the [EntropicCone](https://github.com/blegat/EntropicCone.jl) package.

The `SDDP` algorithm can be run from any node of the lattice of problems using the following function:
```@docs
SDDP(root::SDDPNode, num_stages; K::Int, stopcrit::AbstractStoppingCriterion, verbose, pathsel::Symbol, ztol)
```

This lattice can be built from a [StructJuMP](github.com/StructJuMP/StructJuMP.jl) model using the following function:
```@docs
model2lattice(m::JuMP.Model, num_stages, solver, pruningalgo::CutPruners.AbstractCutPruningAlgo, cutmode::Symbol, newcut::Symbol)
```

## Index

```@index
```
