mutable struct FiniteTempBasis{S, K} <: AbstractBasis{S}
	ptr::Ptr{spir_basis}
	kernel::K
	sve_result::SVEResult{K}
	beta::Float64
	wmax::Float64
	epsilon::Float64
	function FiniteTempBasis{S}(kernel::K, sve_result::SVEResult{K}, β::Real, ωmax::Real, ε::Real) where {S<:Statistics, K<:AbstractKernel}
	    # Create basis
	    status = Ref{Int32}(-100)
	    basis = LibSparseIR.spir_basis_new(_statistics_to_c(S), β, ωmax, kernel.ptr, sve_result.ptr, status)
	    status[] == LibSparseIR.SPIR_COMPUTATION_SUCCESS || error("Failed to create FiniteTempBasis")
	    result = new{S, K}(basis, kernel, sve_result, Float64(β), Float64(ωmax), Float64(ε))
	    finalizer(b -> spir_basis_release(b.ptr), result)
	    return result
	end
end

function FiniteTempBasis{S}(β::Real, ωmax::Real, ε::Real; kernel=LogisticKernel(β * ωmax), sve_result=SVEResult(kernel, ε)) where {S<:Statistics}
	FiniteTempBasis{S}(kernel, sve_result, Float64(β), Float64(ωmax), Float64(ε))
end

# Convenience constructor - matches SparseIR.jl signature
function FiniteTempBasis(stat::S, β::Real, ωmax::Real, ε::Real; kernel=LogisticKernel(β * ωmax), sve_result=SVEResult(kernel, ε)) where {S<:Statistics}
    FiniteTempBasis{typeof(stat)}(β, ωmax, ε; kernel, sve_result)
end

# Basis function type
mutable struct BasisFunction
    ptr::Ptr{spir_funcs}
    basis::FiniteTempBasis  # Keep reference to prevent GC
end
