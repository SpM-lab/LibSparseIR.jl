mutable struct FiniteTempBasis{S} <: AbstractBasis{S}
	ptr::Ptr{spir_basis}
	beta::Float64
	wmax::Float64
	epsilon::Float64
	function FiniteTempBasis{S}(β::Real, ωmax::Real, ε::Real) where {S<:Statistics}
	    # Create logistic kernel
	    kernel_status = Ref{Int32}(-100)
	    kernel = LibSparseIR.spir_logistic_kernel_new(β * ωmax, kernel_status)
	    kernel_status[] == LibSparseIR.SPIR_COMPUTATION_SUCCESS || error("Failed to create Kernel in FiniteTempBasis")

	    # Create SVE result
	    sve_status = Ref{Int32}(-100)
	    sve = LibSparseIR.spir_sve_result_new(kernel, ε, sve_status)
	    sve_status[] == LibSparseIR.SPIR_COMPUTATION_SUCCESS || error("Failed to create SVEResult in FiniteTempBasis")

	    # Create basis
	    status = Ref{Int32}(-100)
	    basis = LibSparseIR.spir_basis_new(_statistics_to_c(S), β, ωmax, kernel, sve, status)
	    status[] == LibSparseIR.SPIR_COMPUTATION_SUCCESS || error("Failed to create FiniteTempBasis")
	    result = new{S}(basis, Float64(β), Float64(ωmax), Float64(ε))
	    finalizer(b -> spir_basis_release(b.ptr), result)
	    return result
	end

	function FiniteTempBasis{S}(kernel::AbstractKernel, β::Real, ωmax::Real, ε::Real) where {S<:Statistics}
	    # Create SVE result
	    sve_status = Ref{Int32}(-100)
	    sve = LibSparseIR.spir_sve_result_new(kernel.ptr, ε, sve_status)
	    sve_status[] == LibSparseIR.SPIR_COMPUTATION_SUCCESS || error("Failed to create SVEResult in FiniteTempBasis")

	    # Create basis
	    status = Ref{Int32}(-100)
	    basis = LibSparseIR.spir_basis_new(_statistics_to_c(S), β, ωmax, kernel.ptr, sve, status)
	    status[] == LibSparseIR.SPIR_COMPUTATION_SUCCESS || error("Failed to create FiniteTempBasis")
	    result = new{S}(basis, Float64(β), Float64(ωmax), Float64(ε))
	    finalizer(b -> spir_basis_release(b.ptr), result)
	    return result
	end
end

# Convenience constructor - matches SparseIR.jl signature
function FiniteTempBasis(stat::S, β::Real, ωmax::Real, ε=nothing; max_size=nothing, kernel=nothing, sve_result=nothing) where {S<:Statistics}
    # Handle optional ε parameter like SparseIR.jl
    if ε === nothing
        # Use default epsilon
        ε = 1e-10  # Same default as in SparseIR.jl
    end

    # Note: max_size, kernel, and sve_result parameters are currently ignored
    # as the C API handles these internally
    # TODO: Add support for these parameters in future versions

    FiniteTempBasis{typeof(stat)}(β, ωmax, ε)
end

# Constructor with explicit kernel parameter - matches SparseIR.jl signature
function FiniteTempBasis(stat::S, kernel::AbstractKernel, β::Real, ωmax::Real, ε=nothing; max_size=nothing, sve_result=nothing) where {S<:Statistics}
    # Validate RegularizedBoseKernel can only be used with Bosonic statistics
    if kernel isa RegularizedBoseKernel && stat isa Fermionic
        throw(ArgumentError("RegularizedBoseKernel does not support fermionic functions"))
    end

    # Handle optional ε parameter like SparseIR.jl
    if ε === nothing
        # Use default epsilon
        ε = 1e-10  # Same default as in SparseIR.jl
    end

    # Note: max_size and sve_result parameters are currently ignored
    # as the C API handles these internally
    # TODO: Add support for these parameters in future versions

    FiniteTempBasis{typeof(stat)}(kernel, β, ωmax, ε)
end

# Basis function type
mutable struct BasisFunction
    ptr::Ptr{spir_funcs}
    basis::FiniteTempBasis  # Keep reference to prevent GC
end
