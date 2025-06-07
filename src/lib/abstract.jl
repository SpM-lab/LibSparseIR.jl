"""
    AbstractBasis

Abstract base class for bases on the imaginary-time axis.

Let `basis` be an abstract basis. Then we can expand a two-point
propagator  `G(τ)`, where `τ` is imaginary time, into a set of basis
functions:

    G(τ) == sum(basis.u[l](τ) * g[l] for l in 1:length(basis)) + ϵ(τ),

where `basis.u[l]` is the `l`-th basis function, `g[l]` is the associated
expansion coefficient and `ϵ(τ)` is an error term. Similarly, the Fourier
transform `Ĝ(n)`, where `n` is now a Matsubara frequency, can be expanded
as follows:

    Ĝ(n) == sum(basis.uhat[l](n) * g[l] for l in 1:length(basis)) + ϵ(n),

where `basis.uhat[l]` is now the Fourier transform of the basis function.
"""
abstract type AbstractBasis{S<:Statistics} end

@doc raw"""
    AbstractKernel

Integral kernel `K(x, y)`.

Abstract base type for an integral kernel, i.e. a AbstractFloat binary function
``K(x, y)`` used in a Fredhold integral equation of the first kind:
```math
    u(x) = ∫ K(x, y) v(y) dy
```
where ``x ∈ [x_\mathrm{min}, x_\mathrm{max}]`` and
``y ∈ [y_\mathrm{min}, y_\mathrm{max}]``. For its SVE to exist,
the kernel must be square-integrable, for its singular values to decay
exponentially, it must be smooth.

In general, the kernel is applied to a scaled spectral function ``ρ'(y)`` as:
```math
    ∫ K(x, y) ρ'(y) dy,
```
where ``ρ'(y) = w(y) ρ(y)``.
"""
abstract type AbstractKernel end

abstract type AbstractReducedKernel <: AbstractKernel end

###############################################################################

"""
    AbstractSVEHints

Discretization hints for singular value expansion of a given kernel.
"""
abstract type AbstractSVEHints end

###############################################################################

"""
    AbstractSampling

Abstract type for sparse sampling.

Encodes the "basis transformation" of a propagator from the truncated IR
basis coefficients `G_ir[l]` to time/frequency sampled on sparse points
`G(x[i])` together with its inverse, a least squares fit:

         ________________                   ___________________
        |                |    evaluate     |                   |
        |     Basis      |---------------->|     Value on      |
        |  coefficients  |<----------------|  sampling points  |
        |________________|      fit        |___________________|
"""
abstract type AbstractSampling{T,Tmat,F} end

###############################################################################

abstract type AbstractSVE end
