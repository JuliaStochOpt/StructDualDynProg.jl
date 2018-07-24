abstract type AbstractStats end

mutable struct Stats <: AbstractStats
    # number of iterations
    niterations::Int
    # n forwards passes of last computation of upper-bound:
    npaths::Int
    # current lower bound
    lowerbound::Float64
    # current lower bound
    upperbound::Float64
    # upper-bound std:
    σ_UB::Float64
end

Stats() = Stats(0, 0, 0.0, Inf, 0.0)

# Add two `Stats` objects and return a new `Stats`
# The second Stats is supposed to be more up to date that the first one
# WARNING: the addition is currently non commutative !
function Base.:+(a::Stats, b::Stats)
    Stats(a.niterations   + b.niterations, b.npaths,
    b.lowerbound, b.upperbound,
    b.σ_UB)
end

function Base.show(io::IO, stats::Stats)
    println(io, "Lower Bound: $(stats.lowerbound)")
    println(io, "Upper Bound: $(stats.upperbound)")
end
