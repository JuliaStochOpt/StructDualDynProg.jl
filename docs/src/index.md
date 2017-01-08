# StochasticDualDynamicProgramming.jl Documentation

## Functions

```@docs
model2lattice(m::JuMP.Model, num_stages, solver, cutmanager::AbstractCutManager, cutmode::Symbol, newcut::Symbol)
SDDP(root::SDDPNode, num_stages; mccount::Int, verbose, pereiracoef, stopcrit::Function, pathsel::Symbol, ztol)
```

## Index

```@index
```
