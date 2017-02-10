# Stochastic Dual Dynamic Programming (SDDP)

| **Documentation** | **Build Status** | **Social** |
|:-----------------:|:----------------:|:----------:|
| | [![Build Status][build-img]][build-url] | [![Gitter][gitter-img]][gitter-url] |
| [![][docs-latest-img]][docs-latest-url] | [![Coveralls branch][coveralls-img]][coveralls-url] [![Codecov branch][codecov-img]][codecov-url] | [<img src="https://upload.wikimedia.org/wikipedia/en/a/af/Discourse_logo.png" width="64">][discourse-url] |

Implementation of [Stochastic Dual Dynamic Programming (SDDP)](http://www.optimization-online.org/DB_FILE/2009/12/2509.pdf).
The problem can either be provided using the [StructJuMP](https://github.com/joehuchette/StructJuMP.jl) modeling interface or using a lower level interface.

This package is used by the [Entropic Cone](https://github.com/blegat/EntropicCone.jl) package.

# Installation
Neither [StructJuMP](https://github.com/joehuchette/StructJuMP.jl) nor this package are currently registered so you will need to use `Pkg.clone` to use them:

```
> Pkg.clone("https://github.com/StructJuMP/StructJuMP.jl")
> Pkg.clone("https://github.com/blegat/StochasticDualDynamicProgramming.jl")
```

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-stable-url]: https://blegat.github.io/StochasticDualDynamicProgramming.jl/stable
[docs-latest-url]: https://blegat.github.io/StochasticDualDynamicProgramming.jl/latest
[build-img]: https://travis-ci.org/blegat/StochasticDualDynamicProgramming.jl.svg?branch=master
[build-url]: https://travis-ci.org/blegat/StochasticDualDynamicProgramming.jl
[coveralls-img]: https://coveralls.io/repos/github/blegat/StochasticDualDynamicProgramming.jl/badge.svg
[coveralls-url]: https://coveralls.io/github/blegat/StochasticDualDynamicProgramming.jl
[codecov-img]: https://codecov.io/gh/blegat/StochasticDualDynamicProgramming.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/blegat/StochasticDualDynamicProgramming.jl
[gitter-url]: https://gitter.im/JuliaOpt/StochasticDualDynamicProgramming.jl
[gitter-img]: https://badges.gitter.im/JuliaOpt/StochasticDualDynamicProgramming.jl.svg
[discourse-url]: https://discourse.julialang.org/c/domain/opt
