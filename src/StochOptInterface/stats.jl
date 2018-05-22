abstract type AbstractSDDPStats end

mutable struct SDDPStats <: AbstractSDDPStats
    # number of calls to solver
    nsolved::Int
    # total time passed inside solver
    solvertime::Float64
    # number of cuts merged
    nmerged::Int
    # time required to merge
    mergetime::Float64
    # number of feasibility cuts
    nfcuts::Int
    fcutstime::Float64
    # number of optimality cuts
    nocuts::Int
    ocutstime::Float64
    nsetx::Int
    setxtime::Float64

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
    # total time
    time::Float64
end

SDDPStats() = SDDPStats(0, 0.0, 0, 0.0, 0, 0.0, 0, 0.0, 0, 0.0, 0, 0, 0, Inf, 0, 0)

# Add two `SDDPStats` objects and return a new `SDDPStats`
# The second SDDPStats is supposed to be more up to date that the first one
# WARNING: the addition is currently non commutative !
function Base.:(+)(a::SDDPStats, b::SDDPStats)
    SDDPStats(a.nsolved + b.nsolved, a.solvertime + b.solvertime,
    a.nmerged + b.nmerged, a.mergetime  + b.mergetime,
    a.nfcuts  + b.nfcuts , a.fcutstime  + b.fcutstime,
    a.nocuts  + b.nocuts , a.ocutstime  + b.ocutstime,
    a.nsetx   + b.nsetx  , a.setxtime   + b.setxtime,
    a.niterations   + b.niterations  , b.npaths,
    b.lowerbound, b.upperbound,
    b.σ_UB, a.time + b.time)
end

macro _time(x)
    quote
        y = @timed $(esc(x))
        # y[1] is returned value
        # y[2] is time in seconds
        y[2]
    end
end


function showtime(t::Float64)
    if !isfinite(t)
        "   min    s    ms    μs"
    else
        s = Int(floor(t))
        t = (t - s) * 1000
        min = div(s, 60)
        s = mod(s, 60)
        ms = Int(floor(t))
        t = (t - ms) * 1000
        μs = Int(floor(t))
        @sprintf "%3dmin %3ds %3dms %3dμs" min s ms μs
    end
end

function Base.show(io::IO, stat::SDDPStats)
    println(io, "                        |     Total time          |  Number  | Average time")
    @printf io "        Solving problem | %s | %8d | %s\n" showtime(stat.solvertime) stat.nsolved showtime(stat.solvertime / stat.nsolved)
    @printf io "          Merging paths | %s | %8d | %s\n" showtime(stat.mergetime ) stat.nmerged showtime(stat.mergetime  / stat.nmerged)
    @printf io "Adding feasibility cuts | %s | %8d | %s\n" showtime(stat.fcutstime ) stat.nfcuts  showtime(stat.fcutstime  / stat.nfcuts )
    @printf io "Adding  optimality cuts | %s | %8d | %s\n" showtime(stat.ocutstime ) stat.nocuts  showtime(stat.ocutstime  / stat.nocuts )
    @printf io "Setting parent solution | %s | %8d | %s\n" showtime(stat.setxtime  ) stat.nsetx   showtime(stat.setxtime   / stat.nsetx  )
    @printf io "                        | %s |" showtime(stat.solvertime+stat.fcutstime+stat.ocutstime+stat.setxtime)
end
