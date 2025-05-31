mutable struct FiniteTempBasis{S} <: AbstractBasis{S}
	ptr::Ptr{spir_basis}
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
	    result = new{S}(basis)
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
	    result = new{S}(basis)
	    finalizer(b -> spir_basis_release(b.ptr), result)
	    return result
	end
end
