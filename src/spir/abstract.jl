_get_ptr(basis::AbstractBasis) = basis.ptr

Base.broadcastable(b::AbstractBasis) = Ref(b)
Base.firstindex(::AbstractBasis) = 1
Base.length(basis::AbstractBasis) = length(basis.s)

"""
    accuracy(basis::AbstractBasis)

Accuracy of the basis.

Upper bound to the relative error of reprensenting a propagator with
the given number of basis functions (number between 0 and 1).
"""
function accuracy end

"""
    significance(basis::AbstractBasis)

Return vector `σ`, where `0 ≤ σ[i] ≤ 1` is the significance level of the `i`-th
basis function. If `ϵ` is the desired accuracy to which to represent a
propagator, then any basis function where `σ[i] < ϵ` can be neglected.

For the IR basis, we simply have that `σ[i] = s[i] / first(s)`.
"""
function significance end

"""
    s(basis::AbstractBasis)

Get the singular values of the basis.
"""
function s end

"""
    u(basis::AbstractBasis)

Get the u basis functions (imaginary time).
"""
function u end

"""
    v(basis::AbstractBasis)

Get the v basis functions (real frequency).
"""
function v end

"""
    uhat(basis::AbstractBasis)

Get the uhat basis functions (Matsubara frequency).
"""
function uhat end

"""
    default_tau_sampling_points(basis::AbstractBasis)

Default sampling points on the imaginary time/x axis.
"""
function default_tau_sampling_points end

"""
    default_matsubara_sampling_points(basis::AbstractBasis; positive_only=false)

Default sampling points on the imaginary frequency axis.

# Arguments

  - `positive_only::Bool`: Only return non-negative frequencies. This is useful if the
    object to be fitted is symmetric in Matsubura frequency, `ĝ(ω) == conj(ĝ(-ω))`,
    or, equivalently, real in imaginary time.
"""
function default_matsubara_sampling_points end

"""
    statistics(basis::AbstractBasis)

Quantum statistic (Statistics instance, Fermionic() or Bosonic()).
"""
statistics(::AbstractBasis{S}) where {S<:Statistics} = S()

"""
    overlap(a::AbstractVector, basis::AbstractBasis, b::AbstractVector)

Compute the overlap ⟨a|S|b⟩ where S is the singular value matrix of the basis.
"""
function overlap(a::AbstractVector, basis::AbstractBasis, b::AbstractVector)
    length(a) == length(basis) || throw(DimensionMismatch("Length of a must match basis size"))
    length(b) == length(basis) || throw(DimensionMismatch("Length of b must match basis size"))

    svals = s(basis)
    return sum(a[i] * svals[i] * b[i] for i in 1:length(basis))
end

"""
    Λ(basis::AbstractBasis)
    lambda(basis::AbstractBasis)

Basis cutoff parameter, `Λ = β * ωmax`, or None if not present
"""
function Λ end
const lambda = Λ

"""
    ωmax(basis::AbstractBasis)
    wmax(basis::AbstractBasis)

Real frequency cutoff or `nothing` if unscaled basis.
"""
function ωmax end
const wmax = ωmax

"""
    β(basis::AbstractBasis)
    beta(basis::AbstractBasis)

Inverse temperature or `nothing` if unscaled basis.
"""
β(basis::AbstractBasis) = basis.β
const beta = β

"""
    iswellconditioned(basis::AbstractBasis)

Returns true if the sampling is expected to be well-conditioned.
"""
iswellconditioned(::AbstractBasis) = true

###############################################################################

Base.broadcastable(kernel::AbstractKernel) = Ref(kernel)

Base.broadcastable(sampling::AbstractSampling) = Ref(sampling)

function LinearAlgebra.cond(sampling::AbstractSampling)
    cond_num = Ref{Float64}(-1.0)
    status = spir_sampling_get_cond_num(sampling.ptr, cond_num)
    status == SPIR_COMPUTATION_SUCCESS || error("Failed to get condition number: $status")
    return cond_num[]
end

"""
    sampling_points(sampling::AbstractSampling)

Return sampling points.
"""
sampling_points(sampling::AbstractSampling) = sampling.sampling_points

"""
    basis(sampling::AbstractSampling)

Return the IR basis associated with `sampling`.
"""
basis(sampling::AbstractSampling) = sampling.basis

function Base.show(io::IO, ::MIME"text/plain", smpl::S) where {S<:AbstractSampling}
    println(io, "$S with sampling points:")
    for p in sampling_points(smpl)[begin:(end - 1)]
        println(io, " $p")
    end
    print(io, " $(last(sampling_points(smpl)))")
end
