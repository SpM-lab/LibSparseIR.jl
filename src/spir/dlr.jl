# AbstractBasis interface implementation

Base.length(dlr::DiscreteLehmannRepresentation) = length(dlr.poles)
Base.size(dlr::DiscreteLehmannRepresentation) = (length(dlr),)

# Pass through to underlying basis
β(dlr::DiscreteLehmannRepresentation) = β(dlr.basis)
ωmax(dlr::DiscreteLehmannRepresentation) = ωmax(dlr.basis)
Λ(dlr::DiscreteLehmannRepresentation) = Λ(dlr.basis)
accuracy(dlr::DiscreteLehmannRepresentation) = accuracy(dlr.basis)

# DLR-specific methods
sampling_points(dlr::DiscreteLehmannRepresentation) = dlr.poles
significance(dlr::DiscreteLehmannRepresentation) = ones(size(dlr))

function default_tau_sampling_points(dlr::DiscreteLehmannRepresentation)
    default_tau_sampling_points(dlr.basis)
end

function default_matsubara_sampling_points(dlr::DiscreteLehmannRepresentation; kwargs...)
    default_matsubara_sampling_points(dlr.basis; kwargs...)
end

# DLR is not as well-conditioned as IR
iswellconditioned(::DiscreteLehmannRepresentation) = false

# Accessor for the underlying basis
basis(dlr::DiscreteLehmannRepresentation) = dlr.basis

"""
    from_IR(dlr::DiscreteLehmannRepresentation, gl::AbstractArray, dims=1)

Transform from IR basis coefficients to DLR coefficients.

# Arguments
- `dlr`: The DLR basis
- `gl`: IR basis coefficients
- `dims`: Dimension along which the basis coefficients are stored

# Returns
DLR coefficients with the same shape as input, but with size `length(dlr)` along dimension `dims`.
"""
function from_IR(dlr::DiscreteLehmannRepresentation, gl::AbstractArray{T,N}, dims=1) where {T,N}
    # Check dimensions
    size(gl, dims) == length(dlr.basis) || throw(DimensionMismatch("Input array has wrong size along dimension $dims"))

    # Prepare output dimensions
    output_dims = collect(size(gl))
    output_dims[dims] = length(dlr)

    # Determine output type
    output_type = T
    output = Array{output_type,N}(undef, output_dims...)

    # Call appropriate C function
    ndim = N
    input_dims = Int32[size(gl)...]
    target_dim = Int32(dims - 1)  # C uses 0-based indexing
    order = C_API.SPIR_ORDER_COLUMN_MAJOR

    if T <: Real
        ret = C_API.spir_ir2dlr_dd(dlr.ptr, order, ndim, input_dims, target_dim, gl, output)
    elseif T <: Complex
        ret = C_API.spir_ir2dlr_zz(dlr.ptr, order, ndim, input_dims, target_dim, gl, output)
    else
        error("Unsupported type: $T")
    end

    ret == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to transform IR to DLR: status=$ret")
    return output
end

"""
    to_IR(dlr::DiscreteLehmannRepresentation, g_dlr::AbstractArray, dims=1)

Transform from DLR coefficients to IR basis coefficients.

# Arguments
- `dlr`: The DLR basis
- `g_dlr`: DLR coefficients
- `dims`: Dimension along which the DLR coefficients are stored

# Returns
IR basis coefficients with the same shape as input, but with size `length(dlr.basis)` along dimension `dims`.
"""
function to_IR(dlr::DiscreteLehmannRepresentation, g_dlr::AbstractArray{T,N}, dims=1) where {T,N}
    # Check dimensions
    size(g_dlr, dims) == length(dlr) || throw(DimensionMismatch("Input array has wrong size along dimension $dims"))

    # Prepare output dimensions
    output_dims = collect(size(g_dlr))
    output_dims[dims] = length(dlr.basis)

    # Determine output type
    output_type = T
    output = Array{output_type,N}(undef, output_dims...)

    # Call appropriate C function
    ndim = N
    input_dims = Int32[size(g_dlr)...]
    target_dim = Int32(dims - 1)  # C uses 0-based indexing
    order = C_API.SPIR_ORDER_COLUMN_MAJOR

    if T <: Real
        ret = C_API.spir_dlr2ir_dd(dlr.ptr, order, ndim, input_dims, target_dim, g_dlr, output)
    elseif T <: Complex
        ret = C_API.spir_dlr2ir_zz(dlr.ptr, order, ndim, input_dims, target_dim, g_dlr, output)
    else
        error("Unsupported type: $T")
    end

    ret == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to transform DLR to IR: status=$ret")
    return output
end

# Pole access functions

"""
    npoles(dlr::DiscreteLehmannRepresentation)

Get the number of poles in the DLR basis.
"""
function npoles(dlr::DiscreteLehmannRepresentation)
    n_poles = Ref{Int32}(-1)
    ret = C_API.spir_dlr_get_npoles(dlr.ptr, n_poles)
    ret == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to get number of poles")
    return Int(n_poles[])
end

"""
    get_poles(dlr::DiscreteLehmannRepresentation)

Get the pole locations for the DLR basis.

Returns a vector of pole locations on the real-frequency axis.
"""
function get_poles(dlr::DiscreteLehmannRepresentation)
    n = npoles(dlr)
    poles = Vector{Float64}(undef, n)
    ret = C_API.spir_dlr_get_poles(dlr.ptr, poles)
    ret == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to get poles")
    return poles
end

# Convenience function for getting default omega sampling points
"""
    default_omega_sampling_points(basis::AbstractBasis)

Get the default real-frequency sampling points for a basis.

These are the extrema of the highest-order basis function on the real-frequency axis,
which provide near-optimal conditioning for the DLR.
"""
function default_omega_sampling_points(basis::AbstractBasis)
    n_points = Ref{Int32}(-1)
    ret = C_API.spir_basis_get_n_default_ws(basis.ptr, n_points)
    ret == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to get number of default omega points")

    points = Vector{Float64}(undef, n_points[])
    ret = C_API.spir_basis_get_default_ws(basis.ptr, points)
    ret == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to get default omega points")

    return points
end