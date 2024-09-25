module ZernikePolynomials

export ZernikeIndex, NM, OSA, Noll, zernike, zernikecoefficients, evaluatezernike, normalization

# Zernike polynomial-based functions
# Some formulas from Wikipedia, some from: "Standards for Reporting the Optical Aberrations of Eyes", Journal of Refractive Surgery Volume 18 September/October 2002
# Reinier Doelman, 23-12-2018

## Represent the different numbering schemes

"""
    ZernikeIndex

An abstract type for specifying the index of a Zernike polynomial. Valid subtypes are [`NM`](@ref), [`Noll`](@ref), and [`OSA`](@ref).
"""
abstract type ZernikeIndex end

Base.convert(::Type{T}, zi::T) where T <: ZernikeIndex = zi
Base.convert(::Type{T}, zi::ZernikeIndex) where T <: ZernikeIndex = T(zi)

"""
    NM(n::Integer, m::Integer)
    NM(noll::Noll)       # construct from Noll index
    NM(osa::OSA)         # construct from OSA index

The integer pair that defines a Zernike polynomial `Zₙᵐ(ρ,θ)`.

See also: [`Noll`](@ref), [`OSA`](@ref)

# Example:
```julia-repl
julia> NM.(OSA.(0:9))
10-element Vector{NM}:
 NM(0, 0)
 NM(1, -1)
 NM(1, 1)
 NM(2, -2)
 NM(2, 0)
 NM(2, 2)
 NM(3, -3)
 NM(3, -1)
 NM(3, 1)
 NM(3, 3)
"""
struct NM <: ZernikeIndex
    n::Int8
    m::Int8

    function NM(n::Integer, m::Integer)
        (n < abs(m) || isodd(n-m)) && throw(ArgumentError("Invalid Zernike index pair (n,m)=($n,$m)."))
        return new(n, m)
    end
end

"""
    Noll(j::Integer)
    Noll(nm::NM)         # construct from NM index
    Noll(osa::OSA)       # construct from OSA index

The Noll single-index number that defines a Zernike polynomial.

See also: [`NM`](@ref), [`OSA`](@ref)
"""
struct Noll <: ZernikeIndex
    j::Int16

    function Noll(j::Integer)
        j >= 1 || throw(ArgumentError("Invalid Noll index $j."))
        return new(j)
    end
end
function Noll(nm::NM)
    m, n = nm.m, nm.n
    p = if mod(n, 4) ∈ (0, 1)
        m > 0 ? 0 : 1
    else
        m ≥ 0 ? 1 : 0
    end
    return Noll(Int((1//2)*n*(n+1) + abs(m) + p))
end

"""
    OSA(j::Integer)
    OSA(nm::NM)
    OSA(noll::Noll)

The OSA/ANSI single-index number that defines a Zernike polynomial.

See also: [`NM`](@ref), [`Noll`](@ref)
"""
struct OSA <: ZernikeIndex
    j::Int16

    function OSA(j::Integer)
        j >= 0 || throw(ArgumentError("Invalid OSA index $j."))
        return new(j)
    end
end

OSA(nm::NM) = ((m, n) = (nm.m, nm.n); OSA(Int((1//2)*(n*(n+2)+m))))


# Other conversions

function NM(osa::OSA)
    j = osa.j
    n = ceil(Int, (-3 + sqrt(9+8j))/2)
    m = 2j-n*(n+2)
    return NM(n, m)
end

function NM(noll::Noll)
    j = noll.j
    n = Int(ceil((-3 + sqrt(1+8j))/2))
    jr = j - Int(n*(n+1)//2)
    if mod(n,4) ∈ (0,1)
        m1 = jr
        m2 = -(jr-1)
        if iseven(n-m1)
            m = m1
        else
            m = m2
        end
    else # mod(n,4) ∈ (2,3)
        m1 = jr-1
        m2 = -(jr)
        if iseven(n-m1)
            m = m1
        else
            m = m2
        end
    end
    return NM(n, m)
end

Noll(osa::OSA) = Noll(NM(osa))
OSA(noll::Noll) = OSA(NM(noll))

## Zernike polynomials

"""
    R(zi::ZernikeIndex)

Obtain the function ρ -> Rₙ^|m|(ρ), where R is the radial polynomial in Zernike polynomials.

# Example:
```julia-repl
julia> rfunc = R(NM(1,1));

julia> rfunc(0.5)
0.5
```
"""
function R(::Type{T}, nm::NM) where T
    m, n = nm.m, nm.n
    p(s) = ((-1)^s * factorial(n-s)) / T(factorial(s) * factorial(Int(0.5 * (n+abs(m)) - s))
                                        * factorial(Int(0.5 * (n-abs(m)) - s)))
    # round brackets to be a generator instead of a Vector
    f(x) = sum(p(s) * x .^ (n-2s) for s in 0:Int((n-abs(m))/2))
    return f

end
R(::Type{T}, zi::ZernikeIndex) where T = R(T, NM(zi))
R(zi::ZernikeIndex) = R(Float64, zi)

"""
    normalization(zi::ZernikeIndex)

Normalization constant of the zernike polynomial.
"""
function normalization(::Type{T}, nm::NM) where T
    m, n = nm.m, nm.n
    δ(x,y) = ifelse(x==y, one(T), zero(T))
    c = sqrt(T(2)*(n+1) / (one(T) + δ(m,0)))
    return c
end
normalization(::Type{T}, zi::ZernikeIndex) where T = normalization(T, NM(zi))
normalization(zi::ZernikeIndex) = normalization(Float64, zi)

"""
    zernike(zi::ZernikeIndex; coord=:polar)

Obtain the function `(ρ,θ) -> Zₙᵐ(ρ,θ)`, where `Zₙᵐ` is the Zernike polynomial corresponding to `zi`. `ρ` is the radius and `θ` the angle.

If coord=:cartesian the function is `(x,y) -> Zₙᵐ(x,y)` in Cartesian coordinates.

# Example:
```julia-repl
julia> zernike(NM(1,1))
julia> zernike(NM(1,1); coord=:polar)
julia> Z = zernike(NM(1,1); coord=:cartesian)  # Z is a polynomial
julia> Z(0.5,0.2)                              # evaluate the polynomial at given (x, y)
```
"""
function zernike(nm::NM; coord=:polar)
    δ(ρ) = ifelse(abs(ρ) ≤ 1, one(eltype(ρ)), zero(eltype(ρ)))

    m, n = nm.m, nm.n
    # use let block to prevent this bug https://github.com/JuliaLang/julia/issues/15276
    # further, we pass the types to normalization
    Z = let
        if m ≥ 0
            (ρ, θ) -> (  normalization(promote_type(eltype(ρ), eltype(θ)), nm)
                       * R(eltype(ρ), nm)(ρ) * cos(m*θ) * δ(ρ))
        else
            (ρ, θ) -> (- normalization(promote_type(eltype(ρ), eltype(θ)), nm)
                       * R(eltype(ρ), nm)(ρ) * sin(m*θ) * δ(ρ))
        end
    end

    if coord == :cartesian
        g(x,y) = (sqrt(x.^2 + y.^2), atan(y,x))
        return (x,y) -> Z(g(x,y)...)
    elseif coord == :polar
        return Z
    else
        throw(ArgumentError("Unrecognized coordinate system $coord"))
    end
end
zernike(zi::ZernikeIndex; kwargs...) = zernike(NM(zi); kwargs...)

"""
    zernikecoefficients(x::AbstractVector, y::AbstractVector, phase::AbstractMatrix, J::AbstractVector{<:ZernikeIndex})
    zernikecoefficients(x::AbstractVector, phase::AbstractMatrix, J::AbstractVector{<:ZernikeIndex})
    zernikecoefficients(phase::AbstractMatrix, J::AbstractVector{<:ZernikeIndex})

Compute the Zernike coefficients (OSA normalization) that in a least-squares sense optimally describe the phase.
The Zernike polynomials used are specified with an index vector `J`.

The default for `x` and `y` is to span pupil coordinates from -1 to 1. If only `x` is supplied, it is also used for `y`.
"""
function zernikecoefficients(X::AbstractVector{<:Real}, Y::AbstractVector{<:Real}, phase::AbstractArray{T,2}, J::AbstractVector{<:ZernikeIndex}) where T
    (length(X), length(Y)) == size(phase) || throw(DimensionMismatch("Size of phase array $(size(phase)) does not match the grid size $(length(X))×$(length(Y))"))

    G = zeros(T, length(X)*length(Y), length(J))
    i = 0
    for j in J
        Z = zernike(j; coord=:cartesian)
        G[:, i+=1] = vec([Z(x, y) for x in X, y in Y])
    end
    return G \ vec(phase)
end
zernikecoefficients(x::AbstractVector, phase::AbstractArray{T,2}, J::AbstractVector{<:ZernikeIndex}) where T =
    zernikecoefficients(x, x, phase, J)

function zernikecoefficients(phase::AbstractArray{T,2}, J::AbstractVector{<:ZernikeIndex}) where T
    s = size(phase)
    X = range(-one(T), stop=one(T), length=s[1])
    Y = range(-one(T), stop=one(T), length=s[2])
    return zernikecoefficients(X, Y, phase, J)
end

"""
    evaluatezernike(x::AbstractVector, y::AbstractVector, J::AbstractVector{<:ZernikeIndex}, coefficients::AbstractVector)
    evaluatezernike(N::Int, J::AbstractVector{<:ZernikeIndex}, coefficients::AbstractVector)

Evaluate the sum of zernike polynomials on a grid specified by `x` and `y`. The coefficient of Zernike polynomial `J[i]` is `coefficients[i]`.

Alternatively, specify the grid size `N`, and the grid will be `N`x`N` with `x` and `y` spanning from -1 to 1.

# Example:
```julia-repl
julia> W = evaluatezernike(64, OSA.[5, 6], [0.3, 4.1])
```
"""
function evaluatezernike(x::AbstractVector, y::AbstractVector, J::AbstractVector{<:ZernikeIndex}, coefficients::AbstractVector{T}) where T
    length(J) == length(coefficients) || throw(ArgumentError("Length of `J` and `coefficients` must match, got $(length(J)) and $(length(coefficients))"))

    out_arr = zeros(T, length(x), length(y))
    for (j, c) in zip(J, coefficients)
        Z = zernike(j; coord=:cartesian)
        out_arr .+= c .* Z.(x, y')
    end

    return out_arr
end
evaluatezernike(x::AbstractVector, J::AbstractVector{<:ZernikeIndex}, coefficients::AbstractVector{T}) where T =
    evaluatezernike(x, x, J, coefficients)

function evaluatezernike(N::Int, J::AbstractVector{<:ZernikeIndex}, coefficients::AbstractVector{T}) where T
    x = range(-one(T), stop=one(T), length=N)
    y = range(-one(T), stop=one(T), length=N)
    return evaluatezernike(x, y, J, coefficients)
end



function evaluatezernike(x, y, zi::ZernikeIndex, coefficients::Real)
    return evaluatezernike(x, y, [zi], [coefficients])
end
function evaluatezernike(n, zi::ZernikeIndex, coefficients::Real)
    return evaluatezernike(n, [zi], [coefficients])
end

end # module
