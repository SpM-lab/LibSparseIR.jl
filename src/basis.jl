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

# Convenience constructor
function FiniteTempBasis(stat::S, β::Real, ωmax::Real, ε::Real) where {S<:Statistics}
    FiniteTempBasis{typeof(stat)}(β, ωmax, ε)
end

# Override the length function from abstract.jl to work with our C API
function Base.length(basis::FiniteTempBasis)
    size_ref = Ref{Int32}(-1)
    ret = C_API.spir_basis_get_size(basis.ptr, size_ref)
    ret == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to get basis size")
    return Int(size_ref[])
end

# Property accessors
β(basis::FiniteTempBasis) = basis.beta
ωmax(basis::FiniteTempBasis) = basis.wmax
Λ(basis::FiniteTempBasis) = basis.beta * basis.wmax

# For now, accuracy is approximated by epsilon
# In reality, it would be computed from the singular values
accuracy(basis::FiniteTempBasis) = basis.epsilon

# Default sampling points
function default_tau_sampling_points(basis::FiniteTempBasis)
    n_points = Ref{Int32}(-1)
    ret = C_API.spir_basis_get_n_default_taus(basis.ptr, n_points)
    ret == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to get number of default tau points")
    
    points = Vector{Float64}(undef, n_points[])
    ret = C_API.spir_basis_get_default_taus(basis.ptr, points)
    ret == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to get default tau points")
    
    return points
end

function default_matsubara_sampling_points(basis::FiniteTempBasis; positive_only=false)
    n_points = Ref{Int32}(-1)
    ret = C_API.spir_basis_get_n_default_matsus(basis.ptr, positive_only, n_points)
    ret == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to get number of default Matsubara points")
    
    points = Vector{Int64}(undef, n_points[])
    ret = C_API.spir_basis_get_default_matsus(basis.ptr, positive_only, points)
    ret == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to get default Matsubara points")
    
    # Convert to MatsubaraFreq objects based on statistics
    if statistics(basis) isa Fermionic
        return [FermionicFreq(n) for n in points]
    else
        return [BosonicFreq(n) for n in points]
    end
end
