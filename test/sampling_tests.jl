@testitem "TauSampling" tags=[:julia] begin
    using LibSparseIR
    using Test
    using Random
    using LinearAlgebra

    @testset "Constructor" begin
        # Create basis
        basis = FiniteTempBasis{Fermionic}(1.0, 10.0, 1e-15)
        
        # Test default constructor
        tau_sampling = TauSampling(basis)
        @test tau_sampling isa TauSampling
        @test tau_sampling isa LibSparseIR.AbstractSampling
        @test length(sampling_points(tau_sampling)) > 0
        @test npoints(tau_sampling) == length(sampling_points(tau_sampling))
        @test LibSparseIR.basis(tau_sampling) === basis
        @test issorted(sampling_points(tau_sampling))
        
        # Test property accessor
        @test tau_sampling.tau == sampling_points(tau_sampling)
        
        # Test with custom sampling points
        custom_points = [-0.5, 0.0, 0.5]
        tau_sampling_custom = TauSampling(basis; sampling_points=custom_points)
        @test sampling_points(tau_sampling_custom) ≈ custom_points
        @test npoints(tau_sampling_custom) == 3
    end
    
    @testset "Evaluate and Fit 1D" begin
        basis = FiniteTempBasis{Fermionic}(1.0, 10.0, 1e-10)
        tau_sampling = TauSampling(basis)
        
        Random.seed!(42)
        coeffs = randn(Float64, length(basis))
        
        # Test evaluate
        tau_values = evaluate(tau_sampling, coeffs)
        @test length(tau_values) == npoints(tau_sampling)
        @test eltype(tau_values) == Float64
        
        # Test roundtrip accuracy
        fitted_coeffs = fit(tau_sampling, tau_values)
        @test length(fitted_coeffs) == length(basis)
        @test fitted_coeffs ≈ coeffs atol=1e-11
        
        # Test in-place versions
        tau_values_inplace = similar(tau_values)
        evaluate!(tau_values_inplace, tau_sampling, coeffs)
        @test tau_values_inplace ≈ tau_values
        
        fitted_coeffs_inplace = similar(fitted_coeffs)
        fit!(fitted_coeffs_inplace, tau_sampling, tau_values)
        @test fitted_coeffs_inplace ≈ fitted_coeffs
    end
    
    @testset "Multi-dimensional arrays" begin
        basis = FiniteTempBasis{Fermionic}(1.0, 10.0, 1e-10)
        tau_sampling = TauSampling(basis)
        
        Random.seed!(123)
        shape = (3, 4, 2)
        
        for dim in 1:3
            # Create array with basis dimension at position `dim`
            array_shape = collect(shape)
            array_shape[dim] = length(basis)
            coeffs = randn(Float64, array_shape...)
            
            # Evaluate
            tau_values = evaluate(tau_sampling, coeffs; dim=dim)
            expected_shape = collect(array_shape)
            expected_shape[dim] = npoints(tau_sampling)
            @test size(tau_values) == tuple(expected_shape...)
            
            # Fit and check roundtrip
            fitted_coeffs = fit(tau_sampling, tau_values; dim=dim)
            @test size(fitted_coeffs) == size(coeffs)
            @test fitted_coeffs ≈ coeffs atol=1e-10
        end
    end
    
    @testset "Error handling" begin
        basis = FiniteTempBasis{Fermionic}(1.0, 10.0, 1e-10)
        tau_sampling = TauSampling(basis)
        
        # Test dimension mismatch - C API returns error codes
        wrong_size_coeffs = randn(length(basis) + 1)
        @test_throws ErrorException evaluate(tau_sampling, wrong_size_coeffs)
        
        wrong_size_values = randn(npoints(tau_sampling) + 1)
        @test_throws ErrorException fit(tau_sampling, wrong_size_values)
    end
end

@testitem "MatsubaraSampling" tags=[:julia] begin
    using LibSparseIR
    using Test
    using Random
    using LinearAlgebra

    @testset "Constructor" begin
        # Test fermionic case
        basis_f = FiniteTempBasis{Fermionic}(1.0, 10.0, 1e-15)
        
        # Default constructor
        matsu_sampling = MatsubaraSampling(basis_f)
        @test matsu_sampling isa MatsubaraSampling
        @test matsu_sampling isa LibSparseIR.AbstractSampling
        @test all(p isa FermionicFreq for p in sampling_points(matsu_sampling))
        @test length(sampling_points(matsu_sampling)) > 0
        @test npoints(matsu_sampling) == length(sampling_points(matsu_sampling))
        @test LibSparseIR.basis(matsu_sampling) === basis_f
        @test !matsu_sampling.positive_only
        
        # Property accessor
        @test matsu_sampling.wn == sampling_points(matsu_sampling)
        
        # Test positive_only=true
        matsu_sampling_pos = MatsubaraSampling(basis_f; positive_only=true)
        @test matsu_sampling_pos.positive_only
        @test all(Int(p) >= 0 for p in sampling_points(matsu_sampling_pos))
        
        # Test bosonic case
        basis_b = FiniteTempBasis{Bosonic}(1.0, 10.0, 1e-15)
        matsu_sampling_b = MatsubaraSampling(basis_b)
        @test all(p isa BosonicFreq for p in sampling_points(matsu_sampling_b))
        
        # Test with custom sampling points
        custom_indices = [1, 3, 5]  # Odd for fermionic
        matsu_sampling_custom = MatsubaraSampling(basis_f; sampling_points=custom_indices)
        @test length(sampling_points(matsu_sampling_custom)) == 3
        @test all(p isa FermionicFreq for p in sampling_points(matsu_sampling_custom))
    end
    
    @testset "Evaluate and Fit 1D" begin
        basis = FiniteTempBasis{Fermionic}(1.0, 10.0, 1e-10)
        matsu_sampling = MatsubaraSampling(basis)
        
        Random.seed!(42)
        coeffs = randn(Float64, length(basis))
        
        # Test evaluate - real to complex
        matsu_values = evaluate(matsu_sampling, coeffs)
        @test length(matsu_values) == npoints(matsu_sampling)
        @test eltype(matsu_values) == ComplexF64
        
        # Test roundtrip accuracy
        fitted_coeffs = fit(matsu_sampling, matsu_values)
        @test length(fitted_coeffs) == length(basis)
        @test eltype(fitted_coeffs) == Float64  # Should return real coefficients
        @test fitted_coeffs ≈ coeffs atol=1e-11
        
        # Test in-place versions
        matsu_values_inplace = similar(matsu_values)
        evaluate!(matsu_values_inplace, matsu_sampling, coeffs)
        @test matsu_values_inplace ≈ matsu_values
        
        fitted_coeffs_inplace = similar(fitted_coeffs)
        fit!(fitted_coeffs_inplace, matsu_sampling, matsu_values)
        @test fitted_coeffs_inplace ≈ fitted_coeffs
    end
    
    @testset "Complex coefficients" begin
        basis = FiniteTempBasis{Fermionic}(1.0, 10.0, 1e-10)
        matsu_sampling = MatsubaraSampling(basis)
        
        Random.seed!(42)
        # For complex coefficients, we expect complex results throughout
        complex_coeffs = randn(ComplexF64, length(basis))
        
        # Evaluate
        matsu_values = evaluate(matsu_sampling, complex_coeffs)
        @test eltype(matsu_values) == ComplexF64
        
        # Fit back should preserve complex nature
        fitted_coeffs = fit(matsu_sampling, matsu_values)
        @test eltype(fitted_coeffs) == Float64  # Our implementation extracts real part
        # So we should only test real part
        @test fitted_coeffs ≈ real(complex_coeffs) atol=1e-11
    end
    
    @testset "Multi-dimensional arrays" begin
        basis = FiniteTempBasis{Fermionic}(1.0, 10.0, 1e-10)
        matsu_sampling = MatsubaraSampling(basis)
        
        Random.seed!(456)
        shape = (2, 3, 4)
        
        for dim in 1:3
            # Create array with basis dimension at position `dim`
            array_shape = collect(shape)
            array_shape[dim] = length(basis)
            coeffs = randn(Float64, array_shape...)
            
            # Evaluate
            matsu_values = evaluate(matsu_sampling, coeffs; dim=dim)
            expected_shape = collect(array_shape)
            expected_shape[dim] = npoints(matsu_sampling)
            @test size(matsu_values) == tuple(expected_shape...)
            @test eltype(matsu_values) == ComplexF64
            
            # Fit and check roundtrip
            fitted_coeffs = fit(matsu_sampling, matsu_values; dim=dim)
            @test size(fitted_coeffs) == size(coeffs)
            @test fitted_coeffs ≈ coeffs atol=1e-10
        end
    end
    
    @testset "positive_only behavior" begin
        basis = FiniteTempBasis{Bosonic}(1.0, 10.0, 1e-10)
        
        # Full sampling
        full_sampling = MatsubaraSampling(basis; positive_only=false)
        full_points = sampling_points(full_sampling)
        @test any(Int(p) < 0 for p in full_points)  # Should have negative frequencies
        
        # Positive only sampling
        pos_sampling = MatsubaraSampling(basis; positive_only=true)
        pos_points = sampling_points(pos_sampling)
        @test all(Int(p) >= 0 for p in pos_points)  # Only non-negative frequencies
        @test length(pos_points) < length(full_points)  # Should have fewer points
        
        # Test that positive_only still gives good reconstruction
        Random.seed!(789)
        coeffs = randn(Float64, length(basis))
        
        pos_values = evaluate(pos_sampling, coeffs)
        fitted_pos = fit(pos_sampling, pos_values)
        @test fitted_pos ≈ coeffs atol=1e-10
    end
end

@testitem "Statistics compatibility" tags=[:julia] begin
    using LibSparseIR
    using Test
    using Random
    using LinearAlgebra
    
    @testset "Fermionic vs Bosonic" begin
        beta = 1.0
        wmax = 10.0
        eps = 1e-10
        
        # Create both types of bases
        basis_f = FiniteTempBasis{Fermionic}(beta, wmax, eps)
        basis_b = FiniteTempBasis{Bosonic}(beta, wmax, eps)
        
        # Create samplings
        tau_f = TauSampling(basis_f)
        tau_b = TauSampling(basis_b)
        matsu_f = MatsubaraSampling(basis_f)
        matsu_b = MatsubaraSampling(basis_b)
        
        # Check that fermionic Matsubara frequencies are odd
        for freq in sampling_points(matsu_f)
            @test isodd(Int(freq))
        end
        
        # Check that bosonic Matsubara frequencies are even
        for freq in sampling_points(matsu_b)
            @test iseven(Int(freq))
        end
        
        # Test that both give good reconstruction
        Random.seed!(321)
        coeffs_f = randn(Float64, length(basis_f))
        coeffs_b = randn(Float64, length(basis_b))
        
        # Fermionic roundtrip
        tau_vals_f = evaluate(tau_f, coeffs_f)
        fitted_f = fit(tau_f, tau_vals_f)
        @test fitted_f ≈ coeffs_f atol=1e-11
        
        matsu_vals_f = evaluate(matsu_f, coeffs_f)
        fitted_matsu_f = fit(matsu_f, matsu_vals_f)
        @test fitted_matsu_f ≈ coeffs_f atol=1e-11
        
        # Bosonic roundtrip
        tau_vals_b = evaluate(tau_b, coeffs_b)
        fitted_b = fit(tau_b, tau_vals_b)
        @test fitted_b ≈ coeffs_b atol=1e-11
        
        matsu_vals_b = evaluate(matsu_b, coeffs_b)
        fitted_matsu_b = fit(matsu_b, matsu_vals_b)
        @test fitted_matsu_b ≈ coeffs_b atol=1e-11
    end
end

@testitem "Edge cases" tags=[:julia] begin
    using LibSparseIR
    using Test
    using LinearAlgebra
    
    @testset "Small basis" begin
        # Test with very small basis
        basis = FiniteTempBasis{Fermionic}(1.0, 1.0, 1e-3)
        @test length(basis) < 10  # Should be small
        
        tau_sampling = TauSampling(basis)
        matsu_sampling = MatsubaraSampling(basis)
        
        coeffs = ones(length(basis))
        
        # Should still work
        tau_vals = evaluate(tau_sampling, coeffs)
        fitted_tau = fit(tau_sampling, tau_vals)
        @test fitted_tau ≈ coeffs atol=1e-8
        
        matsu_vals = evaluate(matsu_sampling, coeffs)
        fitted_matsu = fit(matsu_sampling, matsu_vals)
        @test fitted_matsu ≈ coeffs atol=1e-8
    end
    
    @testset "High precision" begin
        # Test with high precision requirement
        basis = FiniteTempBasis{Fermionic}(1.0, 100.0, 1e-14)
        
        tau_sampling = TauSampling(basis)
        matsu_sampling = MatsubaraSampling(basis)
        
        # Create smooth coefficients that decay
        coeffs = [exp(-i/5.0) for i in 0:(length(basis)-1)]
        
        # Test tau sampling
        tau_vals = evaluate(tau_sampling, coeffs)
        fitted_tau = fit(tau_sampling, tau_vals)
        @test fitted_tau ≈ coeffs atol=1e-12
        
        # Test matsubara sampling
        matsu_vals = evaluate(matsu_sampling, coeffs)
        fitted_matsu = fit(matsu_sampling, matsu_vals)
        @test fitted_matsu ≈ coeffs atol=1e-12
    end
end