
# Property accessors
β(basis::FiniteTempBasis) = basis.beta
ωmax(basis::FiniteTempBasis) = basis.wmax
Λ(basis::FiniteTempBasis) = basis.beta * basis.wmax

# For now, accuracy is approximated by epsilon
# In reality, it would be computed from the singular values
accuracy(basis::FiniteTempBasis) = basis.epsilon

function (f::BasisFunction)(freq::MatsubaraFreq)
    return f(freq.n)
end

function rescale(basis::FiniteTempBasis{S}, new_beta::Real) where {S}
    # Rescale basis to new temperature
    new_lambda = Λ(basis) * new_beta / β(basis)
    kernel = LogisticKernel(new_lambda)
    return FiniteTempBasis{S}(kernel, new_beta, ωmax(basis), accuracy(basis))
end

# Additional utility functions
function significance(basis::FiniteTempBasis)
    svals = s(basis)
    return svals / svals[1]
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
