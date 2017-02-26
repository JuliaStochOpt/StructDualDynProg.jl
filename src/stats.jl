
type SDDPStats
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
end

SDDPStats() = SDDPStats(0,.0,0,.0,0,.0,0,.0,0,.0)

import Base: +, show

function +(a::SDDPStats, b::SDDPStats)
    SDDPStats(a.nsolved + b.nsolved, a.solvertime + b.solvertime,
    a.nmerged + b.nmerged, a.mergetime  + b.mergetime,
    a.nfcuts  + b.nfcuts , a.fcutstime  + b.fcutstime,
    a.nocuts  + b.nocuts , a.ocutstime  + b.ocutstime,
    a.nsetx   + b.nsetx  , a.setxtime   + b.setxtime)
end

macro mytime(x)
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
    println("                        |     Total time [s]      |  Number  | Average time [s]")
    @printf "        Solving problem | %s | %8d | %s\n" showtime(stat.solvertime) stat.nsolved showtime(stat.solvertime / stat.nsolved)
    @printf "          Merging paths | %s | %8d | %s\n" showtime(stat.mergetime ) stat.nmerged showtime(stat.mergetime  / stat.nmerged)
    @printf "Adding feasibility cuts | %s | %8d | %s\n" showtime(stat.fcutstime ) stat.nfcuts  showtime(stat.fcutstime  / stat.nfcuts )
    @printf "Adding  optimality cuts | %s | %8d | %s\n" showtime(stat.ocutstime ) stat.nocuts  showtime(stat.ocutstime  / stat.nocuts )
    @printf "Setting parent solution | %s | %8d | %s\n" showtime(stat.setxtime  ) stat.nsetx   showtime(stat.setxtime   / stat.nsetx  )
    @printf "                        | %s |" showtime(stat.solvertime+stat.fcutstime+stat.ocutstime+stat.setxtime)
end
