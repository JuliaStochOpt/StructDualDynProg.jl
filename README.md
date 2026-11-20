# StructDualDynProg

[![CI][build-img]][build-url] [![Codecov][codecov-img]][codecov-url] [![Documentation][docs-latest-img]][docs-latest-url] [![Gitter][gitter-img]][gitter-url] [Discourse][discourse-url]

Generic interface for the [Stochastic Dual Dynamic Programming (SDDP) algorithm](http://www.optimization-online.org/DB_FILE/2009/12/2509.pdf) as well as a generic implementation.
The problem can either be provided using the [StructJuMP](https://github.com/joehuchette/StructJuMP.jl) modeling interface or using a lower level interface.

This package is used by the [Entropic Cone](https://github.com/blegat/EntropicCone.jl) package.

## Documentation

- [**STABLE**][docs-stable-url] &mdash; **most recently tagged version of the documentation.**
- [**LATEST**][docs-latest-url] &mdash; *in-development version of the documentation.*

## Installation

This branch supports JuMP 1.x. To install the package from a local checkout, run
```julia
] develop .
```

## Solver requirements

Pass a JuMP optimizer factory, such as `GLPK.Optimizer`, to
`SOI.stochasticprogram`. The solver must provide dual solutions for optimal
subproblems, dual certificates for infeasible subproblems, and primal rays
for unbounded subproblems.

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-stable-url]: https://juliastochopt.github.io/StructDualDynProg.jl/stable
[docs-latest-url]: https://juliastochopt.github.io/StructDualDynProg.jl/dev

[build-img]: https://github.com/JuliaStochOpt/StructDualDynProg.jl/actions/workflows/ci.yml/badge.svg
[build-url]: https://github.com/JuliaStochOpt/StructDualDynProg.jl/actions/workflows/ci.yml
[codecov-img]: https://codecov.io/gh/JuliaStochOpt/StructDualDynProg.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/JuliaStochOpt/StructDualDynProg.jl

[gitter-url]: https://gitter.im/JuliaOpt/StochasticDualDynamicProgramming.jl?utm_source=share-link&utm_medium=link&utm_campaign=share-link
[gitter-img]: https://badges.gitter.im/JuliaOpt/StochasticDualDynamicProgramming.jl.svg
[discourse-url]: https://discourse.julialang.org/c/domain/opt
