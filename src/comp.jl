#const threshold = 1e-8

myeqzero{T<:Real}(x::T, ztol) = x == zero(T)
myeqzero{T<:AbstractFloat}(x::T, ztol) = -ztol < x < ztol
myeqzero{T<:Real}(x::Vector{T}, ztol) = myeqzero(sum(abs(x)), ztol)

myeq{T<:Real}(x::T, y::T) = x == y
# I check with zero because isapprox(0, 1e-100) is false...
# but isapprox(1e-100, 2e-100) should be false
myeq{T<:AbstractFloat}(x::T, y::T, ztol) = (x == zero(T) ? myeqzero(y, ztol) : (y == zero(T) ? myeqzero(x, ztol) : isapprox(x, y)))
myeq{T<:Real, S<:Real}(x::T, y::S) = myeq(promote(x, y)...)
myeq{T<:AbstractFloat}(x::Vector{T}, y::Vector{T}, ztol) = (x == zero(x) ? myeqzero(y, ztol) : (y == zero(y) ? myeqzero(x, ztol) : isapprox(x, y)))

mylt{T<:Real}(x::T, y::T, ztol) = x < y
mylt{T<:AbstractFloat}(x::T, y::T, ztol) = x < y && !myeq(x, y, ztol)
mylt{S<:Real,T<:Real}(x::S, y::T, ztol) = mylt(promote(x, y)...)
# mygt{S<:Real, T<:Real}(x::S, y::T) = mylt(y, x)
# myleq{T<:Real}(x::T, y::T) = x <= y
# myleq{T<:AbstractFloat}(x::T, y::T) = x <= y || myeq(x, y)
# myleq{S<:Real,T<:Real}(x::S, y::T) = myleq(promote(x, y)...)
# mygeq{T<:Real}(x::T, y::T) = myleq(y, x)
# mypos{T<:Real}(x::T) = mygt(x, zero(T))
# myneg{T<:Real}(x::T) = mylt(x, zero(T))
# mynonneg{T<:Real}(x::T) = mygeq(x, zero(T))
# mynonpos{T<:Real}(x::T) = myleq(x, zero(T))
