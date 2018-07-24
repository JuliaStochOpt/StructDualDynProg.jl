mutable struct Result
    # n forwards passes of last computation of upper-bound:
    npaths::Int
    # current lower bound
    lowerbound::Float64
    # current lower bound
    upperbound::Float64
    # upper-bound std:
    σ_UB::Float64
end
Result() = Result(0, 0.0, Inf, 0.0)
function Base.show(io::IO, result::Result)
    println(io, "Lower Bound: $(result.lowerbound)")
    println(io, "Upper Bound: $(result.upperbound)")
end

using TimerOutputs
ncalls(to::TimerOutput, key::String) = haskey(to.inner_timers, key) ? TimerOutputs.ncalls(to[key]) : 0
const FCUTS_KEY = "fcuts"
nfcuts(to) = ncalls(to, FCUTS_KEY)
const OCUTS_KEY = "ocuts"
nocuts(to) = ncalls(to, OCUTS_KEY)
iteration_key(it) = "iteration $it"

struct Info
    results::Vector{Result}
    timer::TimerOutput
end
Info() = Info(Result[], TimerOutput())
niterations(info::Info) = length(info.results)
result(info::Info, it::Int) = info.results[it]
last_result(info::Info) = result(info, niterations(info))
timer(info::Info, it::Int) = info.timer[iteration_key(it)]
last_timer(info::Info) = timer(info, niterations(info))
# Time is in ns
total_time_ns(info::Info) = (time_ns() - info.timer.start_data.time)
total_time(info::Info) = total_time_ns(inf) / 1e9
function Base.show(io::IO, info::Info)
    println(last_result(info))
    println(info.timer)
end

function print_iteration_summary(info::Info)
    println("Iteration $(niterations(info))")
    println(last_timer)
end
