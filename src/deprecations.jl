export mn2OSA, mn2Noll, OSA2mn, Noll2mn, OSA2Noll, Noll2OSA

function makeindex(caller, j; index=:OSA)
    Base.depwarn("$caller no longer supports (j; index=:OSA) pairs. Use the `ZernikeIndex` types instead, e.g., `OSA(j)`.", caller)
    if index === :OSA
        return isa(j, Integer) ? OSA(j) : OSA.(j)
    elseif index === :Noll
        return isa(j, Integer) ? Noll(j) : Noll.(j)
    else
        throw(ArgumentError("Unknown Zernike sequential index $index."))
    end
end

function mn2OSA(m::Integer,n::Integer)
    Base.depwarn("mn2OSA(m, n) is deprecated. Use the `ZernikeIndex` types instead, e.g., `OSA(NM(n, m))`. Note that `n` comes before `m`.", :mn2OSA)
    return OSA(NM(n, m)).j
end
function mn2Noll(m::Integer,n::Integer)
    Base.depwarn("mn2Noll(m, n) is deprecated. Use the `ZernikeIndex` types instead, e.g., `Noll(NM(n, m))`. Note that `n` comes before `m`.", :mn2Noll)
    return Noll(NM(n, m)).j
end
function OSA2mn(j::Int)
    Base.depwarn("OSA2mn(j) is deprecated. Use the `ZernikeIndex` types instead, e.g., `NM(OSA(j))`.", :OSA2mn)
    nm = NM(OSA(j))
    return nm.m, nm.n
end
function Noll2mn(j::Int)
    Base.depwarn("Noll2mn(j) is deprecated. Use the `ZernikeIndex` types instead, e.g., `NM(Noll(j))`.", :Noll2mn)
    nm = NM(Noll(j))
    return nm.m, nm.n
end
function OSA2Noll(j::Int)
    Base.depwarn("OSA2Noll(j) is deprecated. Use the `ZernikeIndex` types instead, e.g., `Noll(OSA(j))`.", :OSA2Noll)
    return Noll(OSA(j)).j
end
function Noll2OSA(j::Int)
    Base.depwarn("Noll2OSA(j) is deprecated. Use the `ZernikeIndex` types instead, e.g., `OSA(Noll(j))`.", :Noll2OSA)
    return OSA(Noll(j)).j
end

@deprecate normalization(T::Type, m::Int, n::Int) normalization(T, NM(n, m))
@deprecate normalization(m::Int, n::Int) normalization(NM(n, m))


@deprecate Zernike(m::Int, n::Int; coord=:polar) zernike(NM(n, m); coord=coord)
Zernike(j::Int; index=:OSA, coord=:polar) = zernike(makeindex(:Zernike, j; index=index); coord=coord)
@deprecate Zernike(zi::ZernikeIndex; coord=:polar) zernike(zi; coord=coord)  # this covers the change in capitalization


function Zernikecoefficients(phase::AbstractArray{T,2}, J::Vector{Int}; index=:OSA) where T
    J = makeindex(:Zernikecoefficients, J; index=index)
    return zernikecoefficients(phase, J)
end
function Zernikecoefficients(X::AbstractArray{<: AbstractFloat,1}, phase::AbstractArray{Float64,2}, J::Vector{Int}; index=:OSA)
    J = makeindex(:Zernikecoefficients, J; index=index)
    return zernikecoefficients(X, phase, J)
end
@deprecate Zernikecoefficients(args...; kwargs...) zernikecoefficients(args...; kwargs...)  # cover change in capitalization


function evaluateZernike(N::Int, J::Vector{Int}, coefficients::AbstractArray{T,1}; index=:OSA) where T
    J = makeindex(:evaluateZernike, J; index=index)
    return evaluatezernike(N, J, coefficients)
end
function evaluateZernike(X::AbstractArray{<: AbstractFloat,1}, J::Vector{Int},
    coefficients::Vector{T}; index=:OSA) where T
    J = makeindex(:evaluateZernike, J; index=index)
    return evaluatezernike(X, J, coefficients)
end
function evaluateZernike(n::Int, J::Int, coefficients::T; index=:OSA) where T
    J = makeindex(:evaluateZernike, J; index=index)
    return evaluatezernike(n, J, coefficients)
end
function evaluateZernike(X::AbstractArray{<: AbstractFloat,1}, J::Int, coefficients::T; index=:OSA) where T
    J = makeindex(:evaluateZernike, J; index=index)
    return evaluatezernike(X, J, coefficients)
end

@deprecate evaluateZernike(args...; kwargs...) evaluatezernike(args...; kwargs...)  # cover change in capitalization


# This is used in tests but was not exported
@deprecate R(m::Int, n::Int) R(NM(n, m))
@deprecate R(T::Type, m::Int, n::Int) R(T, NM(n, m))
