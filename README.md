# ZernikePolynomials.jl

This package provides functionality to work with [Zernike polynomials](https://en.wikipedia.org/wiki/Zernike_polynomials).

It features the following:
- functions to convert between common serial indices of Zernike polynomials
- functions to evaluate the polynomials on a grid
- functions to estimate Zernike coefficients in a least squares sense from a user-provided input.

## Conversion to and between (Noll & OSA/ANSI) sequential indices

This package provides conversion utility functions between three common ways of specifying a Zernike polynomial.
Numbering conventions are implemented as subtypes of the abstract type `ZernikeIndex`.

The first is to use `NM(n, m)` to specify the polynomial `Zₙᵐ(ρ,θ)`, where `n` and `m` are integers with `n ≥ abs(m)` (further, `n - abs(m)` must be even), and `(ρ,θ)` are polar coordinates (radius and angle).

The second (sequential) index is the OSA/ANSI standard index (0,1,2,3,...), specified as `OSA(j)`.

The third is Noll's sequential index (1,2,3,...), specified as `Noll(j)`.

Invalid inputs to the `ZernikeIndex` constructors result in an `ArgumentError`, which ensures that you can only construct valid indexes. You can convert between types by calling the constructor or using `convert`:

```julia-repl
julia> using ZernikePolynomials

julia> nm = NM(3, -3)
NM(3, -3)

julia> OSA(nm)
OSA(6)

julia> Noll(nm)
Noll(9)

julia> NM(OSA(6))
NM(3, -3)

julia> [NM(n, m) for n in 0:4, m in -5:5]
ERROR: ArgumentError: Invalid Zernike index pair (n,m)=(0,-5).
...

julia> Noll(0)
ERROR: ArgumentError: Invalid Noll index 0.
...

julia> supertype(NM)
ZernikeIndex
```

## Evaluating Zernike polynomials
The Zernike polynomials `(ρ,θ) -> Zₙᵐ(ρ,θ)` or `(x,y) -> Zₙᵐ(x,y)` can be obtained as a function as follows:

```julia-repl
>> zernike(NM(1,1))
>> zernike(NM(1,1),coord=:polar)
>> Z = zernike(Noll(5),coord=:cartesian) # 5th polynomial by Noll's numbering
>> Z = zernike(NM(1,1),coord=:cartesian)
>> Z(0.5,0.2) # evaluate the function at cartesian coordinate (0.5,0.2)
```

Polar coordinates are used by default. The functions return 0 for ρ > 1.

Zernike polynomials and affine combination thereof can easily be evaluated on a grid of points using the function `evaluatezernike()`.

For example, to evaluate 0.5*Z_2 + 0.3*Z_3 (by OSA numbering) on a 256x256 grid, use the following
```julia-repl
>> evaluatezernike(256, OSA.([2, 3]), [0.5, 0.3])
```

To evaluate it on a (square) grid with predefined coordinates, use for example
```julia-repl
>> x = LinRange(-2,2,256)
>> ϕ = evaluatezernike(x, OSA.([2, 3]), [0.5, 0.3])
>> using Plots
>> heatmap(ϕ)
```
Here the range ``x`` is a range that gives the x and y coordinates.

### Zernike polynomial normalization
The Zernike polynomials are normalized according to Thibos et al. - "Standards for Reporting the Optical Aberrations of Eyes", i.e.
$$
N_n^m = \sqrt{2 (n+1) / (1 + δ(m,0))}
$$
where $δ(m,0) = 1$ for $m = 0$ and 0 otherwise.
These normalization constants can be obtained by `normalization(m,n)`:
```julia-repl
>> [normalization(NM(OSA(i))) for i in 0:5]
```
Note that the definition on the [Wikipedia page on Zernike polynomials](https://en.wikipedia.org/wiki/Zernike_polynomials) is different. Here the polynomials are normalized between [-1,1].
When using or reporting (estimated) Zernike coefficients, it is important to be aware of which normalization has been used.

## Estimating Zernike coefficients in a least squares sense
A common use for Zernike polynomials is to approximate a given 2D input.
Since Zernike polynomials are often used in optics to approximate a 2D phase, this is the term used in the function documentation.

The function `zernikecoefficients()` estimates (in a least-squares sense) the optimal coefficients for a sum of Zernike polynomials. A vector of (sequential) indices should be provided to specify which coefficients should be estimated

```julia-repl
>> ϕ = evaluatezernike(256, OSA.([2, 3]), [0.5, 0.3])
>> zernikecoefficients(ϕ, OSA.([2, 3])) # ≈ [0.5, 0.3]
```

It is also possible to specify on which coordinates the phase is defined:
```julia-repl
>> x = LinRange(-2,2,256)
>> ϕ = evaluatezernike(x, Noll.([2, 3]), [0.5, 0.3])
>> zernikecoefficients(x, ϕ, Noll.([1, 3])) # ≈ [0.0, 0.3]
```
Note that in this last example the coefficient for the Zernike polynomial with Noll index 1 is estimated, but this is not present in ϕ. Therefore the estimated coefficient is 0.

## Plotting example
```julia
using Plots
using ZernikePolynomials

x = LinRange(-1,1,256)
ϕ = evaluatezernike(x, Noll.([2, 3]), [0.5, 0.3])
f(x,y) = (x^2 + y^2) <= 1 ? 1. : NaN # NaNs are not plotted
heatmap( ϕ .* [f(X,Y) for X in x, Y in x])
```
