# Stochastic Dual Dynamic Programming

[![Build Status](https://travis-ci.org/blegat/StochasticDualDynamicProgramming.jl.svg?branch=master)](https://travis-ci.org/blegat/StochasticDualDynamicProgramming.jl)

Implementation of [Stochastic Dual Dynamic Programming](http://www.optimization-online.org/DB_FILE/2009/12/2509.pdf).
The problem can either be provided using the [StructJuMP](https://github.com/joehuchette/StructJuMP.jl) modeling interface or using a lower level interface.

This package is used by the [Entropic Cone](https://github.com/blegat/EntropicCone.jl) package.

# Installation
To use this package you will need Julia v0.5 or later and you currently need [this fork of StructJuMP](https://github.com/blegat/StructJuMP.jl).
You also need to clone it since it is not registered yet.

```
> Pkg.clone("https://github.com/blegat/StructJuMP.jl")
> Pkg.clone("https://github.com/blegat/StochasticDualDynamicProgramming.jl")
```
