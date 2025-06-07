@testitem "TauSampling" tags=[:julia, :spir] begin
    # TODO: Implement me
end

@testitem "MatsubaraSampling" tags=[:julia, :spir] begin

@testset "default_matsubara_sampling_points" begin
    kernel = LogisticKernel(10.0)
    basis = FiniteTempBasis(Fermionic(), β, ωmax, ε; kernel)
    points = LibSparseIR.default_matsubara_sampling_points(basis)
    @test length(points) > 0
end

end
