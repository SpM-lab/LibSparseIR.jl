using .C_API
using .C_API: c_complex, spir_sampling
using LinearAlgebra: SVD, mul!

"""
    TauSampling{T,B} <: AbstractSampling

Sparse sampling in imaginary time using the C API.

Allows transformation between IR basis coefficients and sampling points in imaginary time.
"""
mutable struct TauSampling{T<:Real,B<:AbstractBasis} <: AbstractSampling{T,Float64,Nothing}
    sampling_ptr::Ptr{spir_sampling}
    sampling_points::Vector{T}
    basis::B
    
    function TauSampling{T,B}(sampling_ptr::Ptr{spir_sampling}, sampling_points::Vector{T}, basis::B) where {T<:Real,B<:AbstractBasis}
        obj = new{T,B}(sampling_ptr, sampling_points, basis)
        finalizer(_release_sampling, obj)
        return obj
    end
end

"""
    MatsubaraSampling{T,B} <: AbstractSampling

Sparse sampling in Matsubara frequencies using the C API.

Allows transformation between IR basis coefficients and sampling points in Matsubara frequencies.
"""
mutable struct MatsubaraSampling{T<:MatsubaraFreq,B<:AbstractBasis} <: AbstractSampling{T,ComplexF64,Nothing}
    sampling_ptr::Ptr{spir_sampling}
    sampling_points::Vector{T}
    positive_only::Bool
    basis::B
    
    function MatsubaraSampling{T,B}(sampling_ptr::Ptr{spir_sampling}, sampling_points::Vector{T}, positive_only::Bool, basis::B) where {T<:MatsubaraFreq,B<:AbstractBasis}
        obj = new{T,B}(sampling_ptr, sampling_points, positive_only, basis)
        finalizer(_release_sampling, obj)
        return obj
    end
end

function _release_sampling(sampling::Union{TauSampling,MatsubaraSampling})
    if sampling.sampling_ptr != C_NULL
        C_API.spir_sampling_release(sampling.sampling_ptr)
        sampling.sampling_ptr = C_NULL
    end
end

# Convenience constructors

"""
    TauSampling(basis::AbstractBasis; sampling_points=nothing)

Construct a `TauSampling` object from a basis. If `sampling_points` is not provided,
the default tau sampling points from the basis are used.
"""
function TauSampling(basis::AbstractBasis; sampling_points=nothing)
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
    
    status = Ref{Int32}(-100)
    sampling_ptr = C_API.spir_tau_sampling_new(basis.ptr, length(sampling_points), sampling_points, status)
    status[] == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to create tau sampling: status=$(status[])")
    sampling_ptr != C_NULL || error("Failed to create tau sampling: null pointer returned")
    
    return TauSampling{Float64,typeof(basis)}(sampling_ptr, sampling_points, basis)
end

"""
    MatsubaraSampling(basis::AbstractBasis; positive_only=false, sampling_points=nothing)

Construct a `MatsubaraSampling` object from a basis. If `sampling_points` is not provided,
the default Matsubara sampling points from the basis are used.

If `positive_only=true`, assumes functions are symmetric in Matsubara frequency.
"""
function MatsubaraSampling(basis::AbstractBasis; positive_only=false, sampling_points=nothing)
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
    
    status = Ref{Int32}(-100)
    sampling_ptr = C_API.spir_matsu_sampling_new(basis.ptr, positive_only, length(indices), indices, status)
    status[] == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to create Matsubara sampling: status=$(status[])")
    sampling_ptr != C_NULL || error("Failed to create Matsubara sampling: null pointer returned")
    
    return MatsubaraSampling{eltype(sampling_points),typeof(basis)}(sampling_ptr, sampling_points, positive_only, basis)
end

# Common interface functions

# sampling_points and basis are already defined in abstract.jl

"""
    npoints(sampling::AbstractSampling)

Get the number of sampling points.
"""
function npoints(sampling::Union{TauSampling,MatsubaraSampling})
    n_points = Ref{Int32}(-1)
    ret = C_API.spir_sampling_get_npoints(sampling.sampling_ptr, n_points)
    ret == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to get number of sampling points")
    return Int(n_points[])
end

# Evaluation and fitting functions

"""
    evaluate(sampling::AbstractSampling, coeffs::AbstractArray; dim=1)

Evaluate basis coefficients at the sampling points using the C API.

For multidimensional arrays, `dim` specifies which dimension corresponds to the basis coefficients.
"""
function evaluate(sampling::Union{TauSampling,MatsubaraSampling}, coeffs::AbstractArray{T,N}; dim=1) where {T,N}
    # Determine output dimensions
    output_dims = collect(size(coeffs))
    output_dims[dim] = npoints(sampling)
    
    # Determine output type based on sampling type
    if sampling isa TauSampling
        # For complex input, TauSampling should produce complex output
        output_type = T
        output = Array{output_type,N}(undef, output_dims...)
        evaluate!(output, sampling, coeffs; dim=dim)
    else # MatsubaraSampling
        output_type = T <: Real ? ComplexF64 : promote_type(ComplexF64, T)
        output = Array{output_type,N}(undef, output_dims...)
        evaluate!(output, sampling, coeffs; dim=dim)
    end
    
    return output
end

"""
    evaluate!(output::AbstractArray, sampling::AbstractSampling, coeffs::AbstractArray; dim=1)

In-place version of [`evaluate`](@ref). Write results to the pre-allocated `output` array.
"""
function evaluate!(output::AbstractArray{Tout,N}, sampling::TauSampling, coeffs::AbstractArray{Tin,N}; dim=1) where {Tout,Tin,N}
    # Check dimensions
    expected_dims = collect(size(coeffs))
    expected_dims[dim] = npoints(sampling)
    size(output) == tuple(expected_dims...) || throw(DimensionMismatch("Output array has wrong dimensions"))
    
    # Prepare arguments for C API
    ndim = N
    input_dims = Int32[size(coeffs)...]
    target_dim = Int32(dim - 1)  # C uses 0-based indexing
    order = C_API.SPIR_ORDER_COLUMN_MAJOR
    
    # Call appropriate C function based on input/output types
    if Tin <: Real && Tout <: Real
        ret = C_API.spir_sampling_eval_dd(sampling.sampling_ptr, order, ndim, input_dims, target_dim, coeffs, output)
    elseif Tin <: Complex && Tout <: Complex
        ret = C_API.spir_sampling_eval_zz(sampling.sampling_ptr, order, ndim, input_dims, target_dim, coeffs, output)
    else
        error("Type combination not yet supported for TauSampling: input=$Tin, output=$Tout")
    end
    
    ret == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to evaluate sampling: status=$ret")
    return output
end

function evaluate!(output::AbstractArray{Tout,N}, sampling::MatsubaraSampling, coeffs::AbstractArray{Tin,N}; dim=1) where {Tout,Tin,N}
    # Check dimensions  
    expected_dims = collect(size(coeffs))
    expected_dims[dim] = npoints(sampling)
    size(output) == tuple(expected_dims...) || throw(DimensionMismatch("Output array has wrong dimensions"))
    
    # Prepare arguments for C API
    ndim = N
    input_dims = Int32[size(coeffs)...]
    target_dim = Int32(dim - 1)  # C uses 0-based indexing
    order = C_API.SPIR_ORDER_COLUMN_MAJOR
    
    # Call appropriate C function based on input/output types
    if Tin <: Real && Tout <: Complex
        ret = C_API.spir_sampling_eval_dz(sampling.sampling_ptr, order, ndim, input_dims, target_dim, coeffs, output)
    elseif Tin <: Complex && Tout <: Complex
        ret = C_API.spir_sampling_eval_zz(sampling.sampling_ptr, order, ndim, input_dims, target_dim, coeffs, output)
    else
        error("Type combination not supported for MatsubaraSampling: input=$Tin, output=$Tout")
    end
    
    ret == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to evaluate sampling: status=$ret")
    return output
end

"""
    fit(sampling::AbstractSampling, values::AbstractArray; dim=1)

Fit basis coefficients from values at sampling points using the C API.

For multidimensional arrays, `dim` specifies which dimension corresponds to the sampling points.
"""
function fit(sampling::Union{TauSampling,MatsubaraSampling}, values::AbstractArray{T,N}; dim=1) where {T,N}
    # Determine output dimensions
    output_dims = collect(size(values))
    output_dims[dim] = length(sampling.basis)
    
    # Determine output type - typically real for coefficients 
    if sampling isa TauSampling
        # For complex input, we need complex output
        output_type = T
    else # MatsubaraSampling
        # For Matsubara sampling, we need to be careful about type matching
        # The C API might expect complex output even for real input
        output_type = T <: Complex ? T : ComplexF64
    end
    
    output = Array{output_type,N}(undef, output_dims...)
    fit!(output, sampling, values; dim=dim)
    
    # For MatsubaraSampling, if we want real coefficients, extract real part
    if sampling isa MatsubaraSampling && T <: Complex && output_type <: Complex
        # The fitted coefficients should be real for physical reasons
        # Extract real part and return as real array
        real_output = Array{real(output_type),N}(undef, output_dims...)
        real_output .= real.(output)
        return real_output
    end
    
    return output
end

"""
    fit!(output::AbstractArray, sampling::AbstractSampling, values::AbstractArray; dim=1)

In-place version of [`fit`](@ref). Write results to the pre-allocated `output` array.
"""
function fit!(output::AbstractArray{Tout,N}, sampling::TauSampling, values::AbstractArray{Tin,N}; dim=1) where {Tout,Tin,N}
    # Check dimensions
    expected_dims = collect(size(values))
    expected_dims[dim] = length(sampling.basis)
    size(output) == tuple(expected_dims...) || throw(DimensionMismatch("Output array has wrong dimensions"))
    
    # Prepare arguments for C API
    ndim = N
    input_dims = Int32[size(values)...]
    target_dim = Int32(dim - 1)  # C uses 0-based indexing
    order = C_API.SPIR_ORDER_COLUMN_MAJOR
    
    # Call appropriate C function
    if Tin <: Real && Tout <: Real
        ret = C_API.spir_sampling_fit_dd(sampling.sampling_ptr, order, ndim, input_dims, target_dim, values, output)
    elseif Tin <: Complex && Tout <: Complex
        ret = C_API.spir_sampling_fit_zz(sampling.sampling_ptr, order, ndim, input_dims, target_dim, values, output)
    else
        error("Type combination not yet supported for TauSampling fit: input=$Tin, output=$Tout")
    end
    
    ret == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to fit sampling: status=$ret")
    return output
end

function fit!(output::AbstractArray{Tout,N}, sampling::MatsubaraSampling, values::AbstractArray{Tin,N}; dim=1) where {Tout,Tin,N}
    # Check dimensions
    expected_dims = collect(size(values))
    expected_dims[dim] = length(sampling.basis)
    size(output) == tuple(expected_dims...) || throw(DimensionMismatch("Output array has wrong dimensions"))
    
    # Prepare arguments for C API
    ndim = N
    input_dims = Int32[size(values)...]
    target_dim = Int32(dim - 1)  # C uses 0-based indexing
    order = C_API.SPIR_ORDER_COLUMN_MAJOR
    
    # Call appropriate C function based on input/output types
    if Tin <: Complex && Tout <: Complex
        # Use complex-to-complex API and then extract real part if needed
        ret = C_API.spir_sampling_fit_zz(sampling.sampling_ptr, order, ndim, input_dims, target_dim, values, output)
    elseif Tin <: Complex && Tout <: Real
        # Create temporary complex output, then extract real part
        temp_output = Array{ComplexF64,N}(undef, size(output)...)
        ret = C_API.spir_sampling_fit_zz(sampling.sampling_ptr, order, ndim, input_dims, target_dim, values, temp_output)
        ret == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to fit sampling: status=$ret")
        output .= real.(temp_output)
        return output
    else
        error("Type combination not supported for MatsubaraSampling fit: input=$Tin, output=$Tout")
    end
    
    ret == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to fit sampling: status=$ret")
    return output
end

# Convenience property accessors (similar to SparseIR.jl)
Base.getproperty(s::TauSampling, p::Symbol) = p === :tau ? sampling_points(s) : getfield(s, p)
Base.getproperty(s::MatsubaraSampling, p::Symbol) = p === :wn ? sampling_points(s) : getfield(s, p)