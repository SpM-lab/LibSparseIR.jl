"""
    TauSampling{T,B} <: AbstractSampling

Sparse sampling in imaginary time using the C API.

Allows transformation between IR basis coefficients and sampling points in imaginary time.
"""
mutable struct TauSampling{T<:Real,B<:AbstractBasis} <: AbstractSampling{T,Float64,Nothing}
    ptr::Ptr{spir_sampling}
    sampling_points::Vector{T}
    basis::B

    function TauSampling{T,B}(ptr::Ptr{spir_sampling}, sampling_points::Vector{T}, basis::B) where {T<:Real,B<:AbstractBasis}
        obj = new{T,B}(ptr, sampling_points, basis)
        finalizer(s->spir_sampling_release(s.ptr), obj)
        return obj
    end
end

"""
    MatsubaraSampling{T,B} <: AbstractSampling

Sparse sampling in Matsubara frequencies using the C API.

Allows transformation between IR basis coefficients and sampling points in Matsubara frequencies.
"""
mutable struct MatsubaraSampling{T<:MatsubaraFreq,B<:AbstractBasis} <: AbstractSampling{T,ComplexF64,Nothing}
    ptr::Ptr{spir_sampling}
    sampling_points::Vector{T}
    positive_only::Bool
    basis::B

    function MatsubaraSampling{T,B}(ptr::Ptr{spir_sampling}, sampling_points::Vector{T}, positive_only::Bool, basis::B) where {T<:MatsubaraFreq,B<:AbstractBasis}
        obj = new{T,B}(ptr, sampling_points, positive_only, basis)
        finalizer(s->spir_sampling_release(s.ptr), obj)
        return obj
    end
end

# Convenience constructors

"""
    TauSampling(basis::AbstractBasis; sampling_points=nothing, factorize=true)

Construct a `TauSampling` object from a basis. If `sampling_points` is not provided,
the default tau sampling points from the basis are used.

The `factorize` parameter matches SparseIR.jl interface but is currently ignored
as factorization is handled internally by the C API.
"""
function TauSampling(basis::AbstractBasis; sampling_points=nothing, factorize=true)
    if sampling_points === nothing
        # Get default tau sampling points from basis
        status = Ref{Int32}(-100)
        n_points = Ref{Int32}(-1)

        ret = C_API.spir_basis_get_n_default_taus(basis.ptr, n_points)
        ret == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to get number of default tau points")

        points_array = Vector{Float64}(undef, n_points[])
        ret = C_API.spir_basis_get_default_taus(basis.ptr, points_array)
        ret == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to get default tau points")

        sampling_points = points_array
    else
        sampling_points = collect(Float64, sampling_points)
    end

    # Note: factorize parameter is currently ignored as the C API handles factorization internally
    # TODO: Add support for factorize parameter in future versions

    status = Ref{Int32}(-100)
    ptr = C_API.spir_tau_sampling_new(basis.ptr, length(sampling_points), sampling_points, status)
    status[] == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to create tau sampling: status=$(status[])")
    ptr != C_NULL || error("Failed to create tau sampling: null pointer returned")

    return TauSampling{Float64,typeof(basis)}(ptr, sampling_points, basis)
end

"""
    MatsubaraSampling(basis::AbstractBasis; positive_only=false, sampling_points=nothing, factorize=true)

Construct a `MatsubaraSampling` object from a basis. If `sampling_points` is not provided,
the default Matsubara sampling points from the basis are used.

If `positive_only=true`, assumes functions are symmetric in Matsubara frequency.
The `factorize` parameter matches SparseIR.jl interface but is currently ignored
as factorization is handled internally by the C API.
"""
function MatsubaraSampling(basis::AbstractBasis; positive_only=false, sampling_points=nothing, factorize=true)
    if sampling_points === nothing
        # Get default Matsubara sampling points from basis
        status = Ref{Int32}(-100)
        n_points = Ref{Int32}(-1)

        ret = C_API.spir_basis_get_n_default_matsus(basis.ptr, positive_only, n_points)
        ret == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to get number of default Matsubara points")

        points_array = Vector{Int64}(undef, n_points[])
        ret = C_API.spir_basis_get_default_matsus(basis.ptr, positive_only, points_array)
        ret == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to get default Matsubara points")

        # Convert to MatsubaraFreq objects based on statistics
        if statistics(basis) isa Fermionic
            sampling_points = [FermionicFreq(n) for n in points_array]
        else
            sampling_points = [BosonicFreq(n) for n in points_array]
        end
    else
        # Convert input to appropriate MatsubaraFreq type
        if statistics(basis) isa Fermionic
            sampling_points = [p isa FermionicFreq ? p : FermionicFreq(Int(p)) for p in sampling_points]
        else
            sampling_points = [p isa BosonicFreq ? p : BosonicFreq(Int(p)) for p in sampling_points]
        end
    end

    # Extract indices for C API
    indices = [Int64(Int(p)) for p in sampling_points]

    # Note: factorize parameter is currently ignored as the C API handles factorization internally
    # TODO: Add support for factorize parameter in future versions

    status = Ref{Int32}(-100)
    ptr = C_API.spir_matsu_sampling_new(basis.ptr, positive_only, length(indices), indices, status)
    status[] == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to create Matsubara sampling: status=$(status[])")
    ptr != C_NULL || error("Failed to create Matsubara sampling: null pointer returned")

    return MatsubaraSampling{eltype(sampling_points),typeof(basis)}(ptr, sampling_points, positive_only, basis)
end

