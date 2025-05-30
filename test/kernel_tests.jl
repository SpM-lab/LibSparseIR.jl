@testitem "kernel.jl" tags=[:julia] begin
	using LibSparseIR

	@testset "Logistic kernel" begin
		lam = 42
		kernel = LogisticKernel(lam)
		@test LibSparseIR.Λ(kernel) == lam
	end

	@testset "Regularized Bose kernel" begin
		lam = 42
		kernel = RegularizedBoseKernel(lam)
		@test LibSparseIR.Λ(kernel) == lam
	end

	@testset "Kernel domain" begin
		status = Ref{Int32}(0)
		k = LibSparseIR.spir_logistic_kernel_new(9.0, status)

		xmin = Ref(0.0)
		xmax = Ref(0.0)
		ymin = Ref(0.0)
		ymax = Ref(0.0)
		status = LibSparseIR.spir_kernel_domain(k, xmin, xmax, ymin, ymax)
		@test status == 0
		@test xmin[] == -1.0
		@test xmax[] == 1.0
		@test ymin[] == -1.0
		@test ymax[] == 1.0
	end
end
