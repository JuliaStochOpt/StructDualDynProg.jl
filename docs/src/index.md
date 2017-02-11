# StochasticDualDynamicProgramming.jl Documentation

This packages aims at providing an implementation of SDDP that is both efficient and modular/flexible.
It features the following:

* 

## Functions

```@docs
model2lattice(m::JuMP.Model, num_stages, solver, cutmanager::AbstractCutManager, cutmode::Symbol, newcut::Symbol)
SDDP(root::SDDPNode, num_stages; mccount::Int, verbose, pereiracoef, stopcrit::Function, pathsel::Symbol, ztol)
```

## Index

```@index
```
