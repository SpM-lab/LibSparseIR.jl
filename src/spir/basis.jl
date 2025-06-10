
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
significance(basis::FiniteTempBasis) = basis.s ./ first(basis.s)

function range_to_length(range::UnitRange)
    isone(first(range)) || error("Range must start at 1.")
    return last(range)
end


"""
    finite_temp_bases(β::Real, ωmax::Real, ε=nothing;
                      kernel=LogisticKernel(β * ωmax), sve_result=SVEResult(kernel; ε))

Construct `FiniteTempBasis` objects for fermion and bosons using the same
`LogisticKernel` instance.
"""
function finite_temp_bases(β::Real, ωmax::Real, ε::Real;
        kernel=LogisticKernel(β * ωmax),
        sve_result=SVEResult(kernel, ε))
    basis_f = FiniteTempBasis{Fermionic}(β, ωmax, ε; sve_result, kernel)
    basis_b = FiniteTempBasis{Bosonic}(β, ωmax, ε; sve_result, kernel)
    return basis_f, basis_b
end
