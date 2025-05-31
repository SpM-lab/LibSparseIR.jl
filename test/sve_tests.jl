@testitem "sve.jl" tags=[:julia] begin
	using LibSparseIR

	@testset "sve_result/LogisticKernel" begin
		kernel = LogisticKernel(10)
		sve_result = LibSparseIR.SVEResult(kernel, 1e-10)
		@test true
	end


	@testset "sve_result/RegularizedBoseKernel" begin
		kernel = RegularizedBoseKernel(10)
		sve_result = LibSparseIR.SVEResult(kernel, 1e-10)
		@test true
	end
end
