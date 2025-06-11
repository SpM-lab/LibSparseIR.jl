function (funcs::Ptr{spir_funcs})(x::Real)
    sz = Ref{Int32}(-1)
    spir_funcs_get_size(funcs, sz) == SPIR_COMPUTATION_SUCCESS || error("Failed to get funcs size")
    ret = Vector{Float64}(undef, Int(sz[]))
    spir_funcs_eval(funcs, x, ret) == SPIR_COMPUTATION_SUCCESS || error("Failed to evaluate funcs")
    return ret
end

function (funcs::Ptr{spir_funcs})(x::Vector{Float64})
    hcat(funcs.(x)...)
end

function Base.getindex(funcs::Ptr{spir_funcs}, i::Int)
    status = Ref{Int32}(-100)
    indices = Vector{Int32}(undef, 1)
    indices[1] = i
    ret = spir_funcs_get_slice(funcs, 1, indices, status)
    status[] == SPIR_COMPUTATION_SUCCESS || error("Failed to get basis function u $status[]")
    return ret
end

Base.getindex(funcs::Ptr{spir_funcs}, I) = [funcs[i] for i in I]

function Base.length(funcs::Ptr{spir_funcs})
    sz = Ref{Int32}(-1)
    spir_funcs_get_size(funcs, sz) == SPIR_COMPUTATION_SUCCESS || error("Failed to get funcs size")
    return Int(sz[])
end

Base.firstindex(funcs::Ptr{spir_funcs}) = 1
Base.lastindex(funcs::Ptr{spir_funcs}) = length(funcs)

function roots(poly)
	nroots_ref = Ref{Int32}(-1)
	LibSparseIR.C_API.spir_funcs_get_n_roots(poly, nroots_ref)
	nroots = nroots_ref[]

	out = Vector{Float64}(undef, nroots)
	
	LibSparseIR.C_API.spir_funcs_get_roots(
		poly, out
	)
	return out
end

function overlap_u(poly, f::F) where {F}
	pts = filter(x -> 0 ≤ x ≤ β, roots(poly))
    q, _ = quadgk(
		x -> poly(x) * f(x),
        unique!(sort!(vcat(pts, [0, β])));
        rtol=eps(), order=10, maxevals=10^4)
	q
end

function overlap_v(poly, f::F) where {F}
	pts = filter(x -> -ωmax ≤ x ≤ ωmax, roots(poly))
    q, _ = quadgk(
		x -> poly(x) * f(x),
        unique!(sort!(vcat(pts, [-ωmax, ωmax])));
        rtol=eps(), order=10, maxevals=10^4)
	q
end

