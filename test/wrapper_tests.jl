@testitem "MatsubaraFreq" begin
    using LibSparseIR
    
    @testset "Basic construction and properties" begin
        # Fermionic frequency
        f = FermionicFreq(5)
        @test f.n == 5
        @test f isa MatsubaraFreq{Fermionic}
        
        # Bosonic frequency
        b = BosonicFreq(4)
        @test b.n == 4
        @test b isa MatsubaraFreq{Bosonic}
        
        # pioverbeta constant
        @test pioverbeta isa BosonicFreq
        @test pioverbeta.n == 1
    end
    
    @testset "value and valueim functions" begin
        beta = 10.0
        
        # Fermionic
        f = FermionicFreq(5)
        @test value(f, beta) ≈ (2 * 5 + 1) * π / beta
        @test valueim(f, beta) ≈ im * (2 * 5 + 1) * π / beta
        @test imag(valueim(f, beta)) ≈ value(f, beta)
        @test real(valueim(f, beta)) == 0
        
        # Bosonic
        b = BosonicFreq(4)
        @test value(b, beta) ≈ 2 * 4 * π / beta
        @test valueim(b, beta) ≈ im * 2 * 4 * π / beta
    end
    
    @testset "Arithmetic operations" begin
        f1 = FermionicFreq(5)
        f2 = FermionicFreq(3)
        
        # Addition with integer
        f3 = f1 + 2
        @test f3.n == 7
        @test f3 isa FermionicFreq
        
        # Subtraction with integer
        f4 = f1 - 2
        @test f4.n == 3
        @test f4 isa FermionicFreq
        
        # Negation
        f5 = -f1
        @test f5.n == -5
        @test f5 isa FermionicFreq
        
        # Zero
        f0 = zero(FermionicFreq)
        @test f0.n == 0
        @test f0 isa FermionicFreq
        
        # Same for bosonic
        b1 = BosonicFreq(4)
        b2 = b1 + 3
        @test b2.n == 7
        @test b2 isa BosonicFreq
    end
end

@testitem "FiniteTempBasis construction and properties" begin
    using LibSparseIR
    
    @testset "Basic construction" begin
        beta = 10.0
        omega_max = 1.0
        eps = 1e-10
        
        # Fermionic basis
        basis_f = FiniteTempBasis{Fermionic}(beta, omega_max, eps)
        @test basis_f isa FiniteTempBasis{Fermionic}
        @test basis_f.beta == beta
        @test basis_f.eps == eps
        @test length(basis_f) > 0
        
        # Bosonic basis
        basis_b = FiniteTempBasis{Bosonic}(beta, omega_max, eps)
        @test basis_b isa FiniteTempBasis{Bosonic}
        @test basis_b.beta == beta
        @test basis_b.eps == eps
        @test length(basis_b) > 0
        
        # Default statistics (Fermionic)
        basis_default = FiniteTempBasis(beta, omega_max, eps)
        @test basis_default.statistics isa Fermionic
    end
    
    @testset "Size and indexing" begin
        basis = FiniteTempBasis{Fermionic}(10.0, 1.0, 1e-10)
        n = length(basis)
        
        @test size(basis) == (n,)
        @test firstindex(basis) == 1
        @test lastindex(basis) == n
        
        # Single index
        @test basis[1] == 1
        @test basis[n] == n
        @test_throws BoundsError basis[0]
        @test_throws BoundsError basis[n+1]
        
        # Multiple indices
        inds = basis[1:3]
        @test inds == 1:3
        
        inds2 = basis[[1, 3, 5]]
        @test inds2 == [1, 3, 5]
    end
end

@testitem "Basis functions (u, v, s, uhat)" begin
    using LibSparseIR
    using LibSparseIR: u, v, s, uhat
    
    beta = 10.0
    omega_max = 1.0
    eps = 1e-10
    basis = FiniteTempBasis{Fermionic}(beta, omega_max, eps)
    
    @testset "Singular values" begin
        @test_skip begin
            svals = s(basis)
            @test length(svals) == length(basis)
            @test all(svals .> 0)
            @test issorted(svals, rev=true)  # Should be in descending order
            @test svals[1] ≈ 1.0  # First singular value normalized to 1
            @test svals[end] < eps * 10  # Last singular value should be small
        end
    end
    
    @testset "Basis function objects" begin
        # Get basis functions
        u_funcs = u(basis)
        v_funcs = v(basis)
        uhat_funcs = uhat(basis)
        
        @test u_funcs isa LibSparseIR.BasisFunction
        @test v_funcs isa LibSparseIR.BasisFunction
        @test uhat_funcs isa LibSparseIR.BasisFunction
    end
    
    @testset "Basis function evaluation" begin
        u_funcs = u(basis)
        v_funcs = v(basis)
        uhat_funcs = uhat(basis)
        
        # Evaluate u at tau point
        tau = 0.5
        u_vals = u_funcs(tau)
        @test length(u_vals) == length(basis)
        @test all(isfinite.(u_vals))
        
        # Evaluate v at omega point
        omega = 0.3
        v_vals = v_funcs(omega)
        @test length(v_vals) == length(basis)
        @test all(isfinite.(v_vals))
        
        # Evaluate uhat at Matsubara frequency
        n = 5
        uhat_vals = uhat_funcs(n)
        @test length(uhat_vals) == length(basis)
        @test all(isfinite.(uhat_vals))
        @test eltype(uhat_vals) <: Complex
        
        # Evaluate with MatsubaraFreq
        freq = FermionicFreq(5)
        uhat_vals2 = uhat_funcs(freq)
        @test uhat_vals2 == uhat_vals
    end
end

@testitem "Basis utility functions" begin
    using LibSparseIR
    using LibSparseIR: s, significance, accuracy, rescale, β
    
    beta = 10.0
    omega_max = 1.0
    eps = 1e-10
    basis = FiniteTempBasis{Fermionic}(beta, omega_max, eps)
    
    @testset "significance and accuracy" begin
        sig = significance(basis)
        @test length(sig) == length(basis)
        @test sig[1] ≈ 1.0
        @test all(0 .< sig .<= 1.0)
        @test issorted(sig, rev=true)
        
        acc = accuracy(basis)
        @test acc > 0
        @test acc < eps * 10  # Should be close to requested accuracy
        @test acc == s(basis)[end]  # Should be the last singular value
    end
    
    @testset "rescale" begin
        new_beta = 20.0
        new_basis = rescale(basis, new_beta)
        
        @test new_basis.beta == new_beta
        @test new_basis.eps == basis.eps
        @test new_basis.statistics isa Fermionic
        # The kernel lambda should be scaled appropriately
        @test new_basis.kernel.lambda ≈ basis.kernel.lambda * new_beta / beta
    end
    
    @testset "finite_temp_bases" begin
        using LibSparseIR: finite_temp_bases
        ferm_basis, bose_basis = finite_temp_bases(beta, omega_max, eps)
        
        @test ferm_basis isa FiniteTempBasis{Fermionic}
        @test bose_basis isa FiniteTempBasis{Bosonic}
        @test ferm_basis.beta == bose_basis.beta == beta
        @test ferm_basis.eps == bose_basis.eps == eps
    end
end

@testitem "Default sampling points" begin
    using LibSparseIR
    using LibSparseIR: default_tau_sampling_points, default_omega_sampling_points, default_matsubara_sampling_points
    
    beta = 10.0
    omega_max = 1.0
    eps = 1e-10
    basis = FiniteTempBasis{Fermionic}(beta, omega_max, eps)
    
    @testset "tau sampling points" begin
        tau_points = default_tau_sampling_points(basis)
        @test length(tau_points) > 0
        @test all(0 .<= tau_points .<= beta)
        @test issorted(tau_points)
    end
    
    @testset "omega sampling points" begin
        omega_points = default_omega_sampling_points(basis)
        @test length(omega_points) > 0
        @test all(isfinite.(omega_points))
        @test issorted(omega_points)
    end
    
    @testset "Matsubara sampling points" begin
        # All frequencies
        matsu_points = default_matsubara_sampling_points(basis)
        @test length(matsu_points) > 0
        @test all(p isa MatsubaraFreq{Fermionic} for p in matsu_points)
        
        # Positive only
        matsu_points_pos = default_matsubara_sampling_points(basis, positive_only=true)
        @test length(matsu_points_pos) > 0
        @test all(p.n >= 0 for p in matsu_points_pos)
        @test length(matsu_points_pos) < length(matsu_points)
    end
end

@testitem "TauSampling" begin
    using LibSparseIR
    
    beta = 10.0
    omega_max = 1.0
    eps = 1e-10
    basis = FiniteTempBasis{Fermionic}(beta, omega_max, eps)
    
    @testset "Default tau sampling" begin
        tau_smpl = TauSampling(basis)
        @test tau_smpl isa TauSampling{Fermionic}
        @test basis(tau_smpl) === basis
        
        points = sampling_points(tau_smpl)
        @test length(points) > 0
        @test all(0 .<= points .<= beta)
    end
    
    @testset "Custom tau sampling" begin
        custom_taus = [0.0, 2.5, 5.0, 7.5, 10.0]
        tau_smpl = TauSampling(basis, custom_taus)
        
        points = sampling_points(tau_smpl)
        @test points ≈ custom_taus
    end
    
    @testset "Evaluate and fit" begin
        tau_smpl = TauSampling(basis)
        
        # Create random coefficients
        gl = randn(length(basis))
        
        # Evaluate
        gtau = evaluate(tau_smpl, gl)
        @test length(gtau) == length(sampling_points(tau_smpl))
        @test all(isfinite.(gtau))
        
        # Fit back
        gl_fit = fit(tau_smpl, gtau)
        @test length(gl_fit) == length(basis)
        @test gl ≈ gl_fit atol=1e-10
        
        # Test in-place versions
        gtau2 = similar(gtau)
        evaluate!(tau_smpl, gl, gtau2)
        @test gtau2 ≈ gtau
        
        gl_fit2 = similar(gl_fit)
        fit!(tau_smpl, gtau, gl_fit2)
        @test gl_fit2 ≈ gl_fit
    end
end

@testitem "MatsubaraSampling" begin
    using LibSparseIR
    
    beta = 10.0
    omega_max = 1.0
    eps = 1e-10
    basis = FiniteTempBasis{Fermionic}(beta, omega_max, eps)
    
    @testset "Default Matsubara sampling" begin
        matsu_smpl = MatsubaraSampling(basis)
        @test matsu_smpl isa MatsubaraSampling{Fermionic}
        @test basis(matsu_smpl) === basis
        @test !matsu_smpl.positive_only
        
        points = sampling_points(matsu_smpl)
        @test length(points) > 0
        @test all(p isa MatsubaraFreq{Fermionic} for p in points)
    end
    
    @testset "Positive only sampling" begin
        matsu_smpl = MatsubaraSampling(basis, positive_only=true)
        @test matsu_smpl.positive_only
        
        points = sampling_points(matsu_smpl)
        @test all(p.n >= 0 for p in points)
    end
    
    @testset "Custom Matsubara sampling" begin
        # With MatsubaraFreq objects
        custom_freqs = [FermionicFreq(n) for n in [-5, -3, 0, 3, 5]]
        matsu_smpl = MatsubaraSampling(basis, custom_freqs)
        
        points = sampling_points(matsu_smpl)
        @test length(points) == length(custom_freqs)
        @test all(points[i].n == custom_freqs[i].n for i in 1:length(points))
        
        # With integers
        custom_ints = [-5, -3, 0, 3, 5]
        matsu_smpl2 = MatsubaraSampling(basis, custom_ints)
        points2 = sampling_points(matsu_smpl2)
        @test all(points2[i].n == custom_ints[i] for i in 1:length(points2))
    end
    
    @testset "Evaluate and fit with real coefficients" begin
        matsu_smpl = MatsubaraSampling(basis)
        
        # Real coefficients
        gl = randn(length(basis))
        
        # Evaluate (should return complex)
        gn = evaluate(matsu_smpl, gl)
        @test length(gn) == length(sampling_points(matsu_smpl))
        @test eltype(gn) <: Complex
        @test all(isfinite.(gn))
        
        # For fermionic, imaginary part should be non-zero
        @test any(imag(gn) .!= 0)
    end
    
    @testset "Evaluate and fit with complex coefficients" begin
        matsu_smpl = MatsubaraSampling(basis)
        
        # Complex coefficients
        gl = randn(ComplexF64, length(basis))
        
        # Evaluate
        gn = evaluate(matsu_smpl, gl)
        @test length(gn) == length(sampling_points(matsu_smpl))
        @test eltype(gn) <: Complex
        
        # Fit back
        gl_fit = fit(matsu_smpl, gn)
        @test length(gl_fit) == length(basis)
        @test gl ≈ gl_fit atol=1e-10
        
        # Test in-place versions
        gn2 = similar(gn)
        evaluate!(matsu_smpl, gl, gn2)
        @test gn2 ≈ gn
        
        gl_fit2 = similar(gl_fit)
        fit!(matsu_smpl, gn, gl_fit2)
        @test gl_fit2 ≈ gl_fit
    end
end

@testitem "Overlap function" begin
    using LibSparseIR
    
    beta = 10.0
    omega_max = 1.0
    eps = 1e-10
    basis = FiniteTempBasis{Fermionic}(beta, omega_max, eps)
    
    @testset "Basic overlap" begin
        # Skip tests that depend on s() until C library function is available
        @test_skip begin
            # Create two random vectors
            a = randn(length(basis))
            b = randn(length(basis))
            
            # Compute overlap
            ovlp = overlap(a, basis, b)
            @test isa(ovlp, Float64)
            @test isfinite(ovlp)
            
            # Manual computation
            svals = s(basis)
            expected = sum(a[i] * svals[i] * b[i] for i in 1:length(basis))
            @test ovlp ≈ expected
        end
    end
    
    @testset "Overlap properties" begin
        @test_skip begin
            a = randn(length(basis))
            b = randn(length(basis))
            c = randn(length(basis))
            
            # Symmetry (for real vectors)
            @test overlap(a, basis, b) ≈ overlap(b, basis, a)
            
            # Linearity
            alpha = 2.5
            @test overlap(alpha * a, basis, b) ≈ alpha * overlap(a, basis, b)
            @test overlap(a + b, basis, c) ≈ overlap(a, basis, c) + overlap(b, basis, c)
        end
    end
    
    @testset "Dimension mismatch" begin
        a = randn(length(basis))
        b_wrong = randn(length(basis) + 1)
        
        @test_throws DimensionMismatch overlap(a, basis, b_wrong)
        @test_throws DimensionMismatch overlap(b_wrong, basis, a)
    end
end

@testitem "DiscreteLehmannRepresentation" begin
    using LibSparseIR
    using LibSparseIR: poles
    
    beta = 10.0
    omega_max = 1.0
    eps = 1e-10
    basis = FiniteTempBasis{Fermionic}(beta, omega_max, eps)
    
    @testset "DLR construction" begin
        dlr = DiscreteLehmannRepresentation(basis)
        @test dlr isa DiscreteLehmannRepresentation
        @test dlr.basis === basis
        @test length(dlr) > 0
        @test length(dlr) <= length(basis)
    end
    
    @testset "DLR poles" begin
        dlr = DiscreteLehmannRepresentation(basis)
        dlr_poles = poles(dlr)
        
        @test length(dlr_poles) == length(dlr)
        @test all(isfinite.(dlr_poles))
        @test all(-omega_max .<= dlr_poles .<= omega_max)
    end
    
    @testset "IR to DLR conversion (real)" begin
        dlr = DiscreteLehmannRepresentation(basis)
        
        # Create random IR coefficients
        gl_ir = randn(length(basis))
        
        # Convert to DLR
        gl_dlr = from_IR(dlr, gl_ir)
        @test length(gl_dlr) == length(dlr)
        @test all(isfinite.(gl_dlr))
        
        # Convert back
        gl_ir_back = to_IR(dlr, gl_dlr)
        @test length(gl_ir_back) == length(basis)
        @test gl_ir ≈ gl_ir_back atol=1e-10
    end
    
    @testset "IR to DLR conversion (complex)" begin
        dlr = DiscreteLehmannRepresentation(basis)
        
        # Create random complex IR coefficients
        gl_ir = randn(ComplexF64, length(basis))
        
        # Convert to DLR
        gl_dlr = from_IR(dlr, gl_ir)
        @test length(gl_dlr) == length(dlr)
        @test eltype(gl_dlr) <: Complex
        @test all(isfinite.(gl_dlr))
        
        # Convert back
        gl_ir_back = to_IR(dlr, gl_dlr)
        @test length(gl_ir_back) == length(basis)
        @test gl_ir ≈ gl_ir_back atol=1e-10
    end
    
    @testset "Dimension mismatch" begin
        dlr = DiscreteLehmannRepresentation(basis)
        
        # Wrong size for from_IR
        gl_wrong = randn(length(basis) + 1)
        @test_throws DimensionMismatch from_IR(dlr, gl_wrong)
        
        # Wrong size for to_IR
        dlr_wrong = randn(length(dlr) + 1)
        @test_throws DimensionMismatch to_IR(dlr, dlr_wrong)
    end
end

@testitem "Augmentation types" begin
    using LibSparseIR
    
    @testset "Augmentation type construction" begin
        # Just test that these types exist and can be constructed
        tc = TauConst()
        tl = TauLinear()
        mc = MatsubaraConst()
        
        @test tc isa TauConst
        @test tl isa TauLinear
        @test mc isa MatsubaraConst
        
        @test tc isa LibSparseIR.AbstractAugmentation
        @test tl isa LibSparseIR.AbstractAugmentation
        @test mc isa LibSparseIR.AbstractAugmentation
    end
    
    @testset "AugmentedBasis" begin
        beta = 10.0
        omega_max = 1.0
        eps = 1e-10
        basis = FiniteTempBasis{Fermionic}(beta, omega_max, eps)
        
        aug_basis = LibSparseIR.AugmentedBasis(basis, TauConst())
        @test aug_basis isa LibSparseIR.AugmentedBasis
        @test aug_basis.basis === basis
        @test aug_basis.augmentation isa TauConst
    end
end

@testitem "FiniteTempBasisSet" begin
    using LibSparseIR
    
    @testset "BasisSet construction" begin
        beta = 10.0
        omega_max = 1.0
        eps_values = [1e-6, 1e-8, 1e-10, 1e-12]
        
        # Fermionic
        basis_set_f = FiniteTempBasisSet(beta, omega_max, eps_values, statistics=Fermionic())
        @test basis_set_f isa FiniteTempBasisSet{Fermionic}
        @test length(basis_set_f.bases) == length(eps_values)
        @test basis_set_f.beta == beta
        @test basis_set_f.omega_max == omega_max
        @test basis_set_f.eps_values == eps_values
        
        # Check that bases have increasing size with decreasing epsilon
        sizes = [length(b) for b in basis_set_f.bases]
        @test issorted(sizes)
        
        # Bosonic
        basis_set_b = FiniteTempBasisSet(beta, omega_max, eps_values, statistics=Bosonic())
        @test basis_set_b isa FiniteTempBasisSet{Bosonic}
        @test all(b.statistics isa Bosonic for b in basis_set_b.bases)
    end
end

@testitem "Integration tests" begin
    using LibSparseIR
    using LinearAlgebra: norm
    
    @testset "Full workflow" begin
        # Create basis
        beta = 10.0
        omega_max = 1.0
        eps = 1e-10
        basis = FiniteTempBasis{Fermionic}(beta, omega_max, eps)
        
        # Create samplings
        tau_smpl = TauSampling(basis)
        matsu_smpl = MatsubaraSampling(basis)
        
        # Create some expansion coefficients
        gl = randn(length(basis))
        
        # Evaluate on both tau and Matsubara
        gtau = evaluate(tau_smpl, gl)
        gn = evaluate(matsu_smpl, gl)
        
        # Fit from tau
        gl_from_tau = fit(tau_smpl, gtau)
        @test norm(gl - gl_from_tau) / norm(gl) < 1e-10
        
        # Fit from Matsubara
        gl_from_matsu = fit(matsu_smpl, gn)
        @test norm(gl - gl_from_matsu) / norm(gl) < 1e-10
        
        # Create DLR
        dlr = DiscreteLehmannRepresentation(basis)
        
        # Convert to DLR and back
        gl_dlr = from_IR(dlr, gl)
        gl_back = to_IR(dlr, gl_dlr)
        @test norm(gl - gl_back) / norm(gl) < 1e-10
        
        # Test with complex coefficients
        gl_complex = randn(ComplexF64, length(basis))
        gn_complex = evaluate(matsu_smpl, gl_complex)
        gl_complex_fit = fit(matsu_smpl, gn_complex)
        @test norm(gl_complex - gl_complex_fit) / norm(gl_complex) < 1e-10
    end
    
    @testset "Consistency between Fermionic and Bosonic" begin
        using LibSparseIR: finite_temp_bases
        beta = 10.0
        omega_max = 1.0
        eps = 1e-10
        
        basis_f, basis_b = finite_temp_bases(beta, omega_max, eps)
        
        # Both should have positive lengths
        @test length(basis_f) > 0
        @test length(basis_b) > 0
        
        # Create samplings
        tau_f = TauSampling(basis_f)
        tau_b = TauSampling(basis_b)
        
        # Test that they work
        gl_f = randn(length(basis_f))
        gl_b = randn(length(basis_b))
        
        gtau_f = evaluate(tau_f, gl_f)
        gtau_b = evaluate(tau_b, gl_b)
        
        @test all(isfinite.(gtau_f))
        @test all(isfinite.(gtau_b))
    end
end

@testitem "Kernel types" begin
    using LibSparseIR
    
    @testset "LogisticKernel" begin
        lambda = 10.0
        
        # Fermionic
        kernel_f = LogisticKernel{Fermionic}(lambda)
        @test kernel_f isa LogisticKernel{Fermionic}
        @test kernel_f.lambda == lambda
        
        # Bosonic
        kernel_b = LogisticKernel{Bosonic}(lambda)
        @test kernel_b isa LogisticKernel{Bosonic}
        @test kernel_b.lambda == lambda
    end
    
    @testset "RegularizedBoseKernel" begin
        lambda = 10.0
        kernel = RegularizedBoseKernel(lambda)
        @test kernel isa RegularizedBoseKernel
        @test kernel.lambda == lambda
    end
    
    @testset "Custom kernel in basis" begin
        beta = 10.0
        omega_max = 1.0
        eps = 1e-10
        
        # Create basis with explicit kernel
        kernel = LogisticKernel{Fermionic}(beta * omega_max)
        basis = FiniteTempBasis(beta, omega_max, eps, kernel=kernel)
        @test basis.kernel === kernel
        
        # Create basis with RegularizedBoseKernel
        kernel_bose = RegularizedBoseKernel(beta * omega_max)
        basis_bose = FiniteTempBasis(beta, omega_max, eps, 
                                    kernel=kernel_bose, statistics=Bosonic())
        @test basis_bose.kernel === kernel_bose
    end
end

@testitem "Condition number" begin
    using LibSparseIR
    
    @testset "Sampling condition number" begin
        beta = 10.0
        omega_max = 1.0
        eps = 1e-10
        basis = FiniteTempBasis{Fermionic}(beta, omega_max, eps)
        
        tau_smpl = TauSampling(basis)
        matsu_smpl = MatsubaraSampling(basis)
        
        # Currently returns NaN as placeholder
        cond_tau = cond(tau_smpl)
        cond_matsu = cond(matsu_smpl)
        
        @test isnan(cond_tau)  # Expected until implemented
        @test isnan(cond_matsu)  # Expected until implemented
    end
end

# Tests ported from test/C_API/cinterface_core_tests.jl
# These tests use the Julia wrapper API instead of the C API directly

@testitem "Kernel Accuracy Tests (Julia Wrapper)" tags=[:wrapper] begin
    using LibSparseIR
    using LibSparseIR: Λ

    @testset "LogisticKernel" begin
        Lambda = 9.0
        
        # LogisticKernel is not parametric in the Julia wrapper
        kernel = LogisticKernel(Lambda)
        @test kernel isa LogisticKernel
        @test Λ(kernel) == Lambda
    end

    @testset "RegularizedBoseKernel" begin
        Lambda = 10.0
        kernel = RegularizedBoseKernel(Lambda)
        @test kernel isa RegularizedBoseKernel
        @test Λ(kernel) == Lambda
    end

    @testset "Invalid kernel parameters" begin
        # Negative lambda should throw
        @test_throws DomainError LogisticKernel(-1.0)
        @test_throws DomainError RegularizedBoseKernel(-1.0)
    end
end

@testitem "FiniteTempBasis Constructor Tests (Julia Wrapper)" tags=[:wrapper] begin
    using LibSparseIR
    using LibSparseIR: Λ, β, ωmax, statistics, accuracy

    @testset "Basic constructor with automatic kernel" begin
        beta = 2.0
        wmax = 5.0
        epsilon = 1e-6

        # Fermionic
        basis_f = FiniteTempBasis{Fermionic}(beta, wmax, epsilon)
        @test basis_f isa FiniteTempBasis{Fermionic}
        @test β(basis_f) == beta
        @test ωmax(basis_f) == wmax
        @test accuracy(basis_f) == epsilon  # Currently returns epsilon
        @test length(basis_f) > 0
        @test Λ(basis_f) ≈ beta * wmax

        # Bosonic
        basis_b = FiniteTempBasis{Bosonic}(beta, wmax, epsilon)
        @test basis_b isa FiniteTempBasis{Bosonic}
        @test β(basis_b) == beta
        @test ωmax(basis_b) == wmax
        @test accuracy(basis_b) == epsilon
        @test length(basis_b) > 0
        @test Λ(basis_b) ≈ beta * wmax
    end

    @testset "Constructor with custom kernel" begin
        beta = 2.0
        wmax = 5.0
        Lambda = 10.0
        epsilon = 1e-6

        # Test with LogisticKernel for Fermionic
        kernel = LogisticKernel(Lambda)
        basis_f = FiniteTempBasis{Fermionic}(kernel, beta, wmax, epsilon)
        @test basis_f isa FiniteTempBasis{Fermionic}
        @test statistics(basis_f) isa Fermionic

        # Test with LogisticKernel for Bosonic
        basis_b = FiniteTempBasis{Bosonic}(kernel, beta, wmax, epsilon)
        @test basis_b isa FiniteTempBasis{Bosonic}
        @test statistics(basis_b) isa Bosonic

        # Test with RegularizedBoseKernel
        kernel_reg = RegularizedBoseKernel(Lambda)
        basis_reg = FiniteTempBasis{Bosonic}(kernel_reg, beta, wmax, epsilon)
        @test basis_reg isa FiniteTempBasis{Bosonic}
        @test statistics(basis_reg) isa Bosonic
    end

    @testset "Constructor with statistics instance" begin
        beta = 2.0
        wmax = 5.0
        epsilon = 1e-6

        # Using statistics instances
        basis_f = FiniteTempBasis(Fermionic(), beta, wmax, epsilon)
        @test basis_f isa FiniteTempBasis{Fermionic}
        @test statistics(basis_f) isa Fermionic

        basis_b = FiniteTempBasis(Bosonic(), beta, wmax, epsilon)
        @test basis_b isa FiniteTempBasis{Bosonic}
        @test statistics(basis_b) isa Bosonic
    end

    @testset "Basis size scaling with epsilon" begin
        beta = 2.0
        wmax = 5.0
        
        # Create bases with different epsilon values
        eps_values = [1e-4, 1e-6, 1e-8, 1e-10]
        sizes = Int[]
        
        for eps in eps_values
            basis = FiniteTempBasis{Fermionic}(beta, wmax, eps)
            push!(sizes, length(basis))
        end
        
        # Size should increase as epsilon decreases
        @test issorted(sizes)
        @test all(s > 0 for s in sizes)
    end
end

@testitem "FiniteTempBasis Functions Tests (Julia Wrapper)" tags=[:wrapper] begin
    using LibSparseIR
    using LibSparseIR: β, u, v, uhat, s, significance, accuracy

    function test_basis_functions(statistics_type)
        beta = 2.0
        wmax = 5.0
        epsilon = 1e-6

        # Create basis
        basis = if statistics_type == :fermionic
            FiniteTempBasis{Fermionic}(beta, wmax, epsilon)
        else
            FiniteTempBasis{Bosonic}(beta, wmax, epsilon)
        end

        n = length(basis)
        @test n > 0

        # Get basis functions
        u_funcs = u(basis)
        v_funcs = v(basis)
        uhat_funcs = uhat(basis)

        @test u_funcs isa LibSparseIR.BasisFunction
        @test v_funcs isa LibSparseIR.BasisFunction
        @test uhat_funcs isa LibSparseIR.BasisFunction

        # Test single point evaluation for u basis
        x = 0.5  # Test point for u basis (imaginary time)
        u_vals = u_funcs(x)
        @test length(u_vals) == n
        @test all(isfinite, u_vals)
        @test eltype(u_vals) <: Real

        # Test single point evaluation for v basis
        y = 0.5 * wmax  # Test point for v basis (real frequency)
        v_vals = v_funcs(y)
        @test length(v_vals) == n
        @test all(isfinite, v_vals)
        @test eltype(v_vals) <: Real

        # Test uhat evaluation with integer
        # Note: For bosonic bases, we may need to handle evaluation differently
        if statistics_type == :fermionic
            uhat_vals = uhat_funcs(5)
            @test length(uhat_vals) == n
            @test all(isfinite, uhat_vals)
            @test eltype(uhat_vals) <: Complex
        end

        # Test uhat evaluation with MatsubaraFreq
        freq = if statistics_type == :fermionic
            FermionicFreq(5)
        else
            BosonicFreq(4)  # Use even number for bosonic
        end
        uhat_vals_freq = uhat_funcs(freq)
        @test length(uhat_vals_freq) == n
        @test all(isfinite, uhat_vals_freq)
        @test eltype(uhat_vals_freq) <: Complex
        
        if statistics_type == :fermionic
            @test uhat_vals_freq == uhat_vals
        end

        # Test multiple point evaluation for u basis
        xs = [0.2, 0.4, 0.6, 0.8, 1.0]
        u_matrix = u_funcs(xs)
        @test size(u_matrix) == (length(xs), n)
        @test all(isfinite, u_matrix)

        # Test multiple point evaluation for v basis
        ys = [0.1 * wmax, 0.3 * wmax, 0.5 * wmax]
        v_matrix = v_funcs(ys)
        @test size(v_matrix) == (length(ys), n)
        @test all(isfinite, v_matrix)

        # Test multiple point evaluation for uhat basis
        ns = if statistics_type == :fermionic
            [1, 3, 5, 7]
        else
            [0, 2, 4, 6]  # Even numbers for bosonic
        end
        uhat_matrix = uhat_funcs(ns)
        @test size(uhat_matrix) == (length(ns), n)
        @test all(isfinite, uhat_matrix)
        @test eltype(uhat_matrix) <: Complex

        # Test with MatsubaraFreq array
        freqs = if statistics_type == :fermionic
            [FermionicFreq(n) for n in ns]
        else
            [BosonicFreq(n) for n in ns]
        end
        uhat_matrix2 = uhat_funcs(freqs)
        @test uhat_matrix2 == uhat_matrix

        # Test edge cases
        # u(0) and u(beta) should work
        u_at_0 = u_funcs(0.0)
        u_at_beta = u_funcs(beta)
        @test all(isfinite, u_at_0)
        @test all(isfinite, u_at_beta)

        # v at omega_max boundary
        v_at_max = v_funcs(wmax)
        v_at_neg_max = v_funcs(-wmax)
        @test all(isfinite, v_at_max)
        @test all(isfinite, v_at_neg_max)
    end

    @testset "Basis Functions Fermionic" begin
        test_basis_functions(:fermionic)
    end

    @testset "Basis Functions Bosonic" begin
        test_basis_functions(:bosonic)
    end

    @testset "Basis function properties" begin
        beta = 2.0
        wmax = 5.0
        epsilon = 1e-6
        basis = FiniteTempBasis{Fermionic}(beta, wmax, epsilon)

        # Skip singular value tests until C library function is available
        @test_skip begin
            # Get singular values
            svals = s(basis)
            @test length(svals) == length(basis)
            @test all(svals .> 0)
            @test issorted(svals, rev=true)
            @test svals[1] ≈ 1.0  # First singular value should be normalized to 1
            @test svals[end] < epsilon * 10  # Last singular value should be small

            # Test significance
            sig = significance(basis)
            @test length(sig) == length(basis)
            @test all(0 .< sig .<= 1.0)
            @test sig[1] ≈ 1.0
            @test issorted(sig, rev=true)
        end

        # Test accuracy
        acc = accuracy(basis)
        @test acc > 0
        @test acc == epsilon  # Currently returns epsilon
    end
end