# Common interface functions

# sampling_points and basis are already defined in abstract.jl

"""
    npoints(sampling::AbstractSampling)

Get the number of sampling points.
"""
function npoints(sampling::Union{TauSampling,MatsubaraSampling})
    n_points = Ref{Int32}(-1)
    ret = C_API.spir_sampling_get_npoints(sampling.ptr, n_points)
    ret == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to get number of sampling points")
    return Int(n_points[])
end

# Evaluation and fitting functions

"""
    evaluate(sampling::AbstractSampling, al::AbstractArray; dim=1)

Evaluate basis coefficients at the sampling points using the C API.

For multidimensional arrays, `dim` specifies which dimension corresponds to the basis coefficients.
"""
function evaluate(sampling::Union{TauSampling,MatsubaraSampling}, al::AbstractArray{T,N}; dim=1) where {T,N}
    # Determine output dimensions
    output_dims = collect(size(al))
    output_dims[dim] = npoints(sampling)

    # Determine output type based on sampling type
    if sampling isa TauSampling
        # For complex input, TauSampling should produce complex output
        output_type = T
        output = Array{output_type,N}(undef, output_dims...)
        evaluate!(output, sampling, al; dim=dim)
    else # MatsubaraSampling
        output_type = T <: Real ? ComplexF64 : promote_type(ComplexF64, T)
        output = Array{output_type,N}(undef, output_dims...)
        evaluate!(output, sampling, al; dim=dim)
    end

    return output
end

"""
    evaluate!(output::AbstractArray, sampling::AbstractSampling, al::AbstractArray; dim=1)

In-place version of [`evaluate`](@ref). Write results to the pre-allocated `output` array.
"""
function evaluate!(output::AbstractArray{Tout,N}, sampling::TauSampling, al::AbstractArray{Tin,N}; dim=1) where {Tout,Tin,N}
    # Check dimensions
    expected_dims = collect(size(al))
    expected_dims[dim] = npoints(sampling)
    size(output) == tuple(expected_dims...) || throw(DimensionMismatch("Output array has wrong dimensions"))

    # Prepare arguments for C API
    ndim = N
    input_dims = Int32[size(al)...]
    target_dim = Int32(dim - 1)  # C uses 0-based indexing
    order = C_API.SPIR_ORDER_COLUMN_MAJOR

    # Call appropriate C function based on input/output types
    if Tin <: Real && Tout <: Real
        ret = C_API.spir_sampling_eval_dd(sampling.ptr, order, ndim, input_dims, target_dim, al, output)
    elseif Tin <: Complex && Tout <: Complex
        ret = C_API.spir_sampling_eval_zz(sampling.ptr, order, ndim, input_dims, target_dim, al, output)
    else
        error("Type combination not yet supported for TauSampling: input=$Tin, output=$Tout")
    end

    ret == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to evaluate sampling: status=$ret")
    return output
end

function evaluate!(output::AbstractArray{Tout,N}, sampling::MatsubaraSampling, al::AbstractArray{Tin,N}; dim=1) where {Tout,Tin,N}
    # Check dimensions
    expected_dims = collect(size(al))
    expected_dims[dim] = npoints(sampling)
    size(output) == tuple(expected_dims...) || throw(DimensionMismatch("Output array has wrong dimensions"))

    # Prepare arguments for C API
    ndim = N
    input_dims = Int32[size(al)...]
    target_dim = Int32(dim - 1)  # C uses 0-based indexing
    order = C_API.SPIR_ORDER_COLUMN_MAJOR

    # Call appropriate C function based on input/output types
    if Tin <: Real && Tout <: Complex
        ret = C_API.spir_sampling_eval_dz(sampling.ptr, order, ndim, input_dims, target_dim, al, output)
    elseif Tin <: Complex && Tout <: Complex
        ret = C_API.spir_sampling_eval_zz(sampling.ptr, order, ndim, input_dims, target_dim, al, output)
    else
        error("Type combination not supported for MatsubaraSampling: input=$Tin, output=$Tout")
    end

    ret == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to evaluate sampling: status=$ret")
    return output
end

"""
    fit(sampling::AbstractSampling, al::AbstractArray; dim=1)

Fit basis coefficients from values at sampling points using the C API.

For multidimensional arrays, `dim` specifies which dimension corresponds to the sampling points.
"""
function fit(sampling::Union{TauSampling,MatsubaraSampling}, al::AbstractArray{T,N}; dim=1) where {T,N}
    # Determine output dimensions
    output_dims = collect(size(al))
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
    fit!(output, sampling, al; dim=dim)

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
    fit!(output::AbstractArray, sampling::AbstractSampling, al::AbstractArray; dim=1)

In-place version of [`fit`](@ref). Write results to the pre-allocated `output` array.
"""
function fit!(output::AbstractArray{Tout,N}, sampling::TauSampling, al::AbstractArray{Tin,N}; dim=1) where {Tout,Tin,N}
    # Check dimensions
    expected_dims = collect(size(al))
    expected_dims[dim] = length(sampling.basis)
    size(output) == tuple(expected_dims...) || throw(DimensionMismatch("Output array has wrong dimensions"))

    # Prepare arguments for C API
    ndim = N
    input_dims = Int32[size(al)...]
    target_dim = Int32(dim - 1)  # C uses 0-based indexing
    order = C_API.SPIR_ORDER_COLUMN_MAJOR

    # Call appropriate C function
    if Tin <: Real && Tout <: Real
        ret = C_API.spir_sampling_fit_dd(sampling.ptr, order, ndim, input_dims, target_dim, al, output)
    elseif Tin <: Complex && Tout <: Complex
        ret = C_API.spir_sampling_fit_zz(sampling.ptr, order, ndim, input_dims, target_dim, al, output)
    else
        error("Type combination not yet supported for TauSampling fit: input=$Tin, output=$Tout")
    end

    ret == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to fit sampling: status=$ret")
    return output
end

function fit!(output::AbstractArray{Tout,N}, sampling::MatsubaraSampling, al::AbstractArray{Tin,N}; dim=1) where {Tout,Tin,N}
    # Check dimensions
    expected_dims = collect(size(al))
    expected_dims[dim] = length(sampling.basis)
    size(output) == tuple(expected_dims...) || throw(DimensionMismatch("Output array has wrong dimensions"))

    # Prepare arguments for C API
    ndim = N
    input_dims = Int32[size(al)...]
    target_dim = Int32(dim - 1)  # C uses 0-based indexing
    order = C_API.SPIR_ORDER_COLUMN_MAJOR

    # Call appropriate C function based on input/output types
    if Tin <: Complex && Tout <: Complex
        # Use complex-to-complex API and then extract real part if needed
        ret = C_API.spir_sampling_fit_zz(sampling.ptr, order, ndim, input_dims, target_dim, al, output)
    elseif Tin <: Complex && Tout <: Real
        # Create temporary complex output, then extract real part
        temp_output = Array{ComplexF64,N}(undef, size(output)...)
        ret = C_API.spir_sampling_fit_zz(sampling.ptr, order, ndim, input_dims, target_dim, al, temp_output)
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
