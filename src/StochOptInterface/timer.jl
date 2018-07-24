using TimerOutputs

ncalls(to::TimerOutput, key::String) = haskey(to.inner_timers, key) ? TimerOutputs.ncalls(to[key]) : 0
const FCUTS_KEY = "fcuts"
nfcuts(to) = ncalls(to, FCUTS_KEY)
const OCUTS_KEY = "ocuts"
nocuts(to) = ncalls(to, OCUTS_KEY)
