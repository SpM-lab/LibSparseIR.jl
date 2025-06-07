"""
    DiscreteLehmannRepresentation{S,B} <: AbstractBasis{S}

Discrete Lehmann representation (DLR) with poles selected according to extrema of IR.

This type wraps the C API DLR functionality. The DLR basis is a variant of the IR basis
that uses a "sketching" approach - representing functions as a linear combination of
poles on the real-frequency axis:

    G(iv) == sum(a[i] / (iv - w[i]) for i in 1:npoles)

# Fields
- `ptr::Ptr{spir_basis}`: Pointer to the C DLR object
- `basis::B`: The underlying IR basis
- `poles::Vector{Float64}`: Pole locations on the real-frequency axis
"""
mutable struct DiscreteLehmannRepresentation{S<:Statistics,B<:AbstractBasis{S}} <: AbstractBasis{S}
    ptr::Ptr{spir_basis}
    basis::B
    poles::Vector{Float64}

    function DiscreteLehmannRepresentation{S,B}(ptr::Ptr{spir_basis}, basis::B, poles::Vector{Float64}) where {S<:Statistics,B<:AbstractBasis{S}}
        obj = new{S,B}(ptr, basis, poles)
        finalizer(s->spir_basis_release(s.ptr), obj)
        return obj
    end
end

"""
    DiscreteLehmannRepresentation(basis::AbstractBasis, poles=nothing)

Construct a DLR basis from an IR basis.

If `poles` is not provided, uses the default omega sampling points from the IR basis.
"""
function DiscreteLehmannRepresentation(basis::AbstractBasis, poles=nothing)
    if poles === nothing
        # Get default omega sampling points from basis
        n_poles = Ref{Int32}(-1)
        ret = C_API.spir_basis_get_n_default_ws(basis.ptr, n_poles)
        ret == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to get number of default omega points")

        poles = Vector{Float64}(undef, n_poles[])
        ret = C_API.spir_basis_get_default_ws(basis.ptr, poles)
        ret == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to get default omega points")

        # Create DLR using default constructor
        status = Ref{Int32}(-100)
        dlr_ptr = C_API.spir_dlr_new(basis.ptr, status)
        status[] == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to create DLR: status=$(status[])")
        dlr_ptr != C_NULL || error("Failed to create DLR: null pointer returned")
    else
        # Create DLR with user-specified poles
        poles = collect(Float64, poles)
        status = Ref{Int32}(-100)
        dlr_ptr = C_API.spir_dlr_new_with_poles(basis.ptr, length(poles), poles, status)
        status[] == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to create DLR with poles: status=$(status[])")
        dlr_ptr != C_NULL || error("Failed to create DLR with poles: null pointer returned")
    end

    # Get the actual poles from the C API as they might differ
    n_poles_actual = Ref{Int32}(-1)
    ret = C_API.spir_dlr_get_npoles(dlr_ptr, n_poles_actual)
    ret == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to get number of DLR poles")

    poles_actual = Vector{Float64}(undef, n_poles_actual[])
    ret = C_API.spir_dlr_get_poles(dlr_ptr, poles_actual)
    ret == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to get DLR poles")

    return DiscreteLehmannRepresentation{typeof(statistics(basis)),typeof(basis)}(dlr_ptr, basis, poles_actual)
end

