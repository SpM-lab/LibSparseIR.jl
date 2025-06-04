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

function default_omega_sampling_points(basis::FiniteTempBasis)
    n_points = Ref{Int32}(-1)
    ret = C_API.spir_basis_get_n_default_ws(basis.ptr, n_points)
    ret == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to get number of default omega points")
    
    points = Vector{Float64}(undef, n_points[])
    ret = C_API.spir_basis_get_default_ws(basis.ptr, points)
    ret == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to get default omega points")
    
    return points
end

# Basis function type
mutable struct BasisFunction
    ptr::Ptr{spir_funcs}
    basis::FiniteTempBasis  # Keep reference to prevent GC
end

# Evaluate basis function at single point
function (f::BasisFunction)(x::Real)
    out = Vector{Float64}(undef, length(f.basis))
    status = C_API.spir_funcs_eval(f.ptr, Float64(x), out)
    status == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to evaluate basis functions")
    return out
end

# Evaluate basis function at Matsubara frequency (for uhat)
function (f::BasisFunction)(n::Integer)
    out = Vector{ComplexF64}(undef, length(f.basis))
    status = C_API.spir_funcs_eval_matsu(f.ptr, Int64(n), out)
    status == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to evaluate basis functions at Matsubara frequency")
    return out
end

function (f::BasisFunction)(freq::MatsubaraFreq)
    return f(freq.n)
end

# Evaluate at multiple points
function (f::BasisFunction)(xs::AbstractVector{<:Real})
    n_points = length(xs)
    n_funcs = length(f.basis)
    out = Matrix{Float64}(undef, n_points, n_funcs)
    
    # Convert to column-major format expected by C API
    status = C_API.spir_funcs_batch_eval(f.ptr, C_API.SPIR_ORDER_COLUMN_MAJOR, 
                                         n_points, collect(Float64, xs), vec(out))
    status == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to evaluate basis functions at multiple points")
    return out
end

# Evaluate at multiple Matsubara frequencies
function (f::BasisFunction)(ns::AbstractVector{<:Integer})
    n_points = length(ns)
    n_funcs = length(f.basis)
    out = Matrix{ComplexF64}(undef, n_points, n_funcs)
    
    # Convert to column-major format expected by C API
    status = C_API.spir_funcs_batch_eval_matsu(f.ptr, C_API.SPIR_ORDER_COLUMN_MAJOR,
                                               n_points, collect(Int64, ns), vec(out))
    status == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to evaluate basis functions at multiple Matsubara frequencies")
    return out
end

function (f::BasisFunction)(freqs::AbstractVector{<:MatsubaraFreq})
    ns = [freq.n for freq in freqs]
    return f(ns)
end

# Get basis functions u, v, s, uhat
function u(basis::FiniteTempBasis)
    status = Ref{Int32}(-100)
    ptr = C_API.spir_basis_get_u(basis.ptr, status)
    status[] == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to get u basis functions")
    funcs = BasisFunction(ptr, basis)
    finalizer(f -> C_API.spir_funcs_release(f.ptr), funcs)
    return funcs
end

function v(basis::FiniteTempBasis)
    status = Ref{Int32}(-100)
    ptr = C_API.spir_basis_get_v(basis.ptr, status)
    status[] == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to get v basis functions")
    funcs = BasisFunction(ptr, basis)
    finalizer(f -> C_API.spir_funcs_release(f.ptr), funcs)
    return funcs
end

function s(basis::FiniteTempBasis)
    n = length(basis)
    svals = Vector{Float64}(undef, n)
    status = C_API.spir_basis_get_singular_values(basis.ptr, svals)
    status == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to get singular values")
    return svals
end

function uhat(basis::FiniteTempBasis)
    status = Ref{Int32}(-100)
    ptr = C_API.spir_basis_get_uhat(basis.ptr, status)
    status[] == C_API.SPIR_COMPUTATION_SUCCESS || error("Failed to get uhat basis functions")
    funcs = BasisFunction(ptr, basis)
    finalizer(f -> C_API.spir_funcs_release(f.ptr), funcs)
    return funcs
end

# Additional utility functions
function significance(basis::FiniteTempBasis)
    svals = s(basis)
    return svals / svals[1]
end

function rescale(basis::FiniteTempBasis{S}, new_beta::Real) where {S}
    # Rescale basis to new temperature
    new_lambda = Λ(basis) * new_beta / β(basis)
    kernel = LogisticKernel(new_lambda)
    return FiniteTempBasis{S}(kernel, new_beta, ωmax(basis), accuracy(basis))
end

function finite_temp_bases(β::Real, ωmax::Real, ε=nothing; kernel=nothing, sve_result=nothing)
    # Handle optional ε parameter like SparseIR.jl
    if ε === nothing
        ε = 1e-10  # Default epsilon value
    end
    
    # Note: kernel and sve_result parameters are currently ignored
    # as the C API handles these internally
    # TODO: Add support for these parameters in future versions
    
    ferm_basis = FiniteTempBasis{Fermionic}(β, ωmax, ε)
    bose_basis = FiniteTempBasis{Bosonic}(β, ωmax, ε)
    return (ferm_basis, bose_basis)
end
