@testitem "basis.jl" tags=[:julia, :lib] begin
	using LibSparseIR

	β = 2.0
	ωmax = 5.0
	ε = 1e-6
	Λ = β * ωmax
	@testset "FiniteTempBasis{S} for S=$(S)" for S in [Fermionic, Bosonic]
		basis = FiniteTempBasis(S(), β, ωmax, ε)
		@test true
	end

	@testset "FiniteTempBasis{S}/LogisticKernel for S=$(S)" for S in [Fermionic, Bosonic]
		kernel = LogisticKernel(Λ)
		basis = FiniteTempBasis(S(), kernel, β, ωmax, ε)
		@test true
	end

	@testset "FiniteTempBasis{S}/RegularizedBoseKernel for S=$(S)" for S in [Fermionic, Bosonic]
		kernel = RegularizedBoseKernel(Λ)
		if S() isa Fermionic
			@test_throws "RegularizedBoseKernel does not support fermionic functions" FiniteTempBasis(S(), kernel, β, ωmax, ε)
		else
			@test true
		end
	end
end
