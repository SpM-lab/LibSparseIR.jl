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

