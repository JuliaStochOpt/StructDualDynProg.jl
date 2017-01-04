# Stochastic Dual Dynamic Programming (SDDP)

[![Build Status](https://travis-ci.org/blegat/StochasticDualDynamicProgramming.jl.svg?branch=master)](https://travis-ci.org/blegat/StochasticDualDynamicProgramming.jl)
[![Coverage Status](https://coveralls.io/repos/github/blegat/StochasticDualDynamicProgramming.jl/badge.svg)](https://coveralls.io/github/blegat/StochasticDualDynamicProgramming.jl)
[![codecov](https://codecov.io/gh/blegat/StochasticDualDynamicProgramming.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/blegat/StochasticDualDynamicProgramming.jl)

Implementation of [Stochastic Dual Dynamic Programming (SDDP)](http://www.optimization-online.org/DB_FILE/2009/12/2509.pdf).
The problem can either be provided using the [StructJuMP](https://github.com/joehuchette/StructJuMP.jl) modeling interface or using a lower level interface.

This package is used by the [Entropic Cone](https://github.com/blegat/EntropicCone.jl) package.

# Installation
Neither [StructJuMP](https://github.com/joehuchette/StructJuMP.jl) nor this package are currently registered so you will need to use `Pkg.clone` to use them:

```
> Pkg.clone("https://github.com/StructJuMP/StructJuMP.jl")
> Pkg.clone("https://github.com/blegat/StochasticDualDynamicProgramming.jl")
```
