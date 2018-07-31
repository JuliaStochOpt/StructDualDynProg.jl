# StructDualDynProg

| **Documentation** | **PackageEvaluator** | **Build Status** | **Social** |
|:-----------------:|:--------------------:|:----------------:|:----------:|
| [![][docs-stable-img]][docs-stable-url] | [![][pkg-0.6-img]][pkg-0.6-url] | [![Build Status][build-img]][build-url] [![Build Status][winbuild-img]][winbuild-url] | [![Gitter][gitter-img]][gitter-url] |
| [![][docs-latest-img]][docs-latest-url] | [![][pkg-0.7-img]][pkg-0.7-url] | [![Coveralls branch][coveralls-img]][coveralls-url] [![Codecov branch][codecov-img]][codecov-url] | [<img src="https://upload.wikimedia.org/wikipedia/commons/thumb/a/af/Discourse_logo.png/799px-Discourse_logo.png" width="64">][discourse-url] |

Generic interface for the [Stochastic Dual Dynamic Programming (SDDP) algorithm](http://www.optimization-online.org/DB_FILE/2009/12/2509.pdf) as well as a generic implementation.
The problem can either be provided using the [StructJuMP](https://github.com/joehuchette/StructJuMP.jl) modeling interface or using a lower level interface.

This package is used by the [Entropic Cone](https://github.com/blegat/EntropicCone.jl) package.

## Documentation

- [**STABLE**][docs-stable-url] &mdash; **most recently tagged version of the documentation.**
- [**LATEST**][docs-latest-url] &mdash; *in-development version of the documentation.*

## Installation

```julia
> Pkg.update()
> Pkg.add("StructDualDynProg")
```

### Development version

The development version uses [StochOptInterface (SOI)](https://github.com/JuliaStochOpt/StochOptInterface.jl)
which is not released yet so to use the development version of this package,
install it as explained above and then do the following:
```julia
> Pkg.clone("https://github.com/JuliaStochOpt/StochOptInterface.jl.git")
> Pkg.checkout("StructDualDynProg")
```

## Notes for choice of solvers
This package should work with any linear programming (LP) solver [supported by MathProgBase](http://www.juliaopt.org/).
If some subproblems are infeasible, an infeasibility ray will be asked to the solver.
In this case, it is advised to turn presolve off in `CPLEX` (i.e. `CPLEX.CplexSolver(CPX_PARAM_REDUCE=0)`) and to avoid using Clp since it [often cannot find the infeasibility ray](https://projects.coin-or.org/Clp/ticket/79).
If some subproblems are unbounded, an unbounded ray and a feasible solution will be asked for the solver. Again, [avoid using Clp in that case](https://projects.coin-or.org/Clp/ticket/82). Hopefully subproblems shouldn't be unbounded if `detectlb` is left at `true` when calling `model2lattice`.

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-stable-url]: https://juliastochopt.github.io/StructDualDynProg.jl/stable
[docs-latest-url]: https://juliastochopt.github.io/StructDualDynProg.jl/latest

[pkg-0.6-img]: http://pkg.julialang.org/badges/StructDualDynProg_0.6.svg
[pkg-0.6-url]: http://pkg.julialang.org/?pkg=StructDualDynProg
[pkg-0.7-img]: http://pkg.julialang.org/badges/StructDualDynProg_0.7.svg
[pkg-0.7-url]: http://pkg.julialang.org/?pkg=StructDualDynProg

[build-img]: https://travis-ci.org/JuliaStochOpt/StructDualDynProg.jl.svg?branch=master
[build-url]: https://travis-ci.org/JuliaStochOpt/StructDualDynProg.jl
[winbuild-img]: https://ci.appveyor.com/api/projects/status/3lh9e6tujojgodar/branch/master?svg=true
[winbuild-url]: https://ci.appveyor.com/project/blegat/structdualdynprog-jl/branch/master
[coveralls-img]: https://coveralls.io/repos/github/blegat/StructDualDynProg.jl/badge.svg
[coveralls-url]: https://coveralls.io/github/JuliaStochOpt/StructDualDynProg.jl
[codecov-img]: https://codecov.io/gh/JuliaStochOpt/StructDualDynProg.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/JuliaStochOpt/StructDualDynProg.jl

[gitter-url]: https://gitter.im/JuliaOpt/StochasticDualDynamicProgramming.jl?utm_source=share-link&utm_medium=link&utm_campaign=share-link
[gitter-img]: https://badges.gitter.im/JuliaOpt/StochasticDualDynamicProgramming.jl.svg
[discourse-url]: https://discourse.julialang.org/c/domain/opt
