# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     custom_cell_magics: kql
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: Julia 1.11.5
#     language: julia
#     name: julia-1.11
# ---

# %%
using Pkg
Pkg.activate(joinpath(@__DIR__, "LibSparseIR.jl"))
using LibSparseIR
using Test

# %%
# helper function corresponding to _get_dims in cinterface_integration.cxx
function _get_dims(target_dim_size::Integer, extra_dims::Vector{<:Integer}, target_dim::Integer, ndim::Integer)
    dims = Vector{Int32}(undef, ndim)
    dims[target_dim+1] = target_dim_size  # Julia is 1-indexed
    pos = 1
    for i in 1:ndim
        if i == target_dim + 1
            continue
        end
        dims[i] = extra_dims[pos]
        pos += 1
    end
    return dims
end

# %%
function generate_random_coeffs(::Type{<:Real}, random_value, pole)
    (2 * random_value - 1.0) * sqrt(abs(pole))
end

function generate_random_coeffs(::Type{<:Complex}, random_value, pole)
    (2 * random_value - 1.0) * sqrt(abs(pole)) + im * (2 * random_value - 1.0) * sqrt(abs(pole))
end

# %%
# Helper function corresponding to _spir_basis_new in cinterface_integration.cxx
function _spir_basis_new(statistics::Integer, beta::Float64, omega_max::Float64, epsilon::Float64, status::Ref{Cint})
    # Create logistic kernel
    kernel_status = Ref{Int32}(0)
    kernel = LibSparseIR.spir_logistic_kernel_new(beta * omega_max, kernel_status)
    @test kernel_status[] == LibSparseIR.SPIR_COMPUTATION_SUCCESS
    @test kernel != C_NULL

    # Create SVE result
    sve_status = Ref{Int32}(0)
    sve = LibSparseIR.spir_sve_result_new(kernel, epsilon, sve_status)
    @test sve_status[] == LibSparseIR.SPIR_COMPUTATION_SUCCESS
    @test sve != C_NULL

    # Create basis
    basis = LibSparseIR.spir_basis_new(statistics, beta, omega_max, kernel, sve, status)
    @test status[] == LibSparseIR.SPIR_COMPUTATION_SUCCESS
    @test basis != C_NULL

    # Clean up intermediate objects (like C++ version)
    LibSparseIR.spir_sve_result_release(sve)
    LibSparseIR.spir_kernel_release(kernel)

    return basis
end

# %%
function getperm(N, src, dst)
    perm = collect(1:N)
    deleteat!(perm, src)
    insert!(perm, dst, src)
    return perm
end

"""
    movedim(arr::AbstractArray, src => dst)

Move `arr`'s dimension at `src` to `dst` while keeping the order of the remaining
dimensions unchanged.
"""
function movedim(arr::AbstractArray{T,N}, src, dst) where {T,N}
    src == dst && return arr
    return permutedims(arr, getperm(N, src, dst))
end

coeffs = movedim(coeffs_targetdim0, 1, 1 + target_dim) # Julia is 1-indexed
nothing

# %%
function dlr_to_IR(dlr, order, ndim, dims, target_dim, coeffs::AbstractArray{<:Real}, g_IR::AbstractArray{<:Real})
    LibSparseIR.spir_dlr2ir_dd(dlr, order, ndim, dims, target_dim, coeffs, g_IR)
end

function dlr_to_IR(dlr, order, ndim, dims, target_dim, coeffs::AbstractArray{<:Complex}, g_IR::AbstractArray{<:Complex})
    LibSparseIR.spir_dlr2ir_zz(dlr, order, ndim, dims, target_dim, coeffs, g_IR)
end

function dlr_from_IR(dlr, order, ndim, dims, target_dim, g_IR::AbstractArray{<:Real}, g_DLR_reconst::AbstractArray{<:Real})
    LibSparseIR.spir_ir2dlr_dd(dlr, order, ndim, dims, target_dim, g_IR, g_DLR_reconst)
end

function dlr_from_IR(dlr, order, ndim, dims, target_dim, g_IR::AbstractArray{<:Complex}, g_DLR_reconst::AbstractArray{<:Complex})
    LibSparseIR.spir_ir2dlr_zz(dlr, order, ndim, dims, target_dim, g_IR, g_DLR_reconst)
end

# %%
beta = 1e+4
beta = 10.0
wmax = 2.0
epsilon = 1e-10
tol = 10 * epsilon
extra_dims = [2, 3, 4]
target_dim = 0 # Julia is 1-indexed C/C++ is 0-indexed
ndim = length(extra_dims) + 1
T = Float64
order = LibSparseIR.SPIR_ORDER_COLUMN_MAJOR
positive_only = false

# %%
status = Ref{Cint}(-100)
stat = LibSparseIR.SPIR_STATISTICS_BOSONIC
basis = _spir_basis_new(1, beta, wmax, epsilon, status)
@test status[] == LibSparseIR.SPIR_COMPUTATION_SUCCESS

basis_size_ref = Ref{Cint}(-100)
status = LibSparseIR.spir_basis_get_size(basis, basis_size_ref)
@test status[] == LibSparseIR.SPIR_COMPUTATION_SUCCESS
basis_size = basis_size_ref[]

# %%
status = Ref{Cint}(-100)
dlr = LibSparseIR.spir_dlr_new(basis, status)
@test status[] == LibSparseIR.SPIR_COMPUTATION_SUCCESS

# %%
npoles_ref = Ref{Cint}(-100)
status = LibSparseIR.spir_dlr_get_npoles(dlr, npoles_ref)
@test status[] == LibSparseIR.SPIR_COMPUTATION_SUCCESS
npoles = npoles_ref[]
@test npoles >= basis_size

poles = Vector{Float64}(undef, npoles)
status = LibSparseIR.spir_dlr_get_poles(dlr, poles)
@test status[] == LibSparseIR.SPIR_COMPUTATION_SUCCESS

# %%
extra_size = prod(extra_dims)

# %%
coeffs_targetdim0 = Array{T,4}(undef, npoles, extra_dims...)
coeffs_2d = reshape(coeffs_targetdim0, Int64(npoles), Int64(extra_size))
for i in axes(coeffs_2d, 1)
    for j in axes(coeffs_2d, 2)
        coeffs_2d[i, j] = generate_random_coeffs(T, rand(), poles[i])
    end
end

# %%
@test poles .|> abs |> maximum <= wmax

# %%
coeffs = movedim(coeffs_targetdim0, 1, 1 + target_dim) # Julia is 1-indexed
nothing

# %%
# Convert DLR coefficients to IR coefficients
g_IR = Array{T,4}(undef, _get_dims(npoles, extra_dims, target_dim, ndim)...)
status = dlr_to_IR(dlr, order, ndim, _get_dims(npoles, extra_dims, target_dim, ndim), target_dim, coeffs, g_IR)
@test status == LibSparseIR.SPIR_COMPUTATION_SUCCESS

# %%
# Convert IR coefficients to DLR coefficients
g_DLR_reconst = Array{T,4}(undef, _get_dims(npoles, extra_dims, target_dim, ndim)...)
status = dlr_from_IR(dlr, order, ndim, _get_dims(npoles, extra_dims, target_dim, ndim), target_dim, g_IR, g_DLR_reconst)
@test status == LibSparseIR.SPIR_COMPUTATION_SUCCESS

# %%
# From IR C API
g_dlr = Array{T,4}(undef, _get_dims(npoles, extra_dims, target_dim, ndim)...)
status = dlr_from_IR(dlr, order, ndim, _get_dims(basis_size, extra_dims, target_dim, ndim), target_dim, g_IR, g_dlr)
@test status == LibSparseIR.SPIR_COMPUTATION_SUCCESS

# %%
# DLR basis functions
dlr_u_status = Ref{Cint}(-100)
dlr_u = LibSparseIR.spir_basis_get_u(dlr, dlr_u_status)
@test dlr_u_status[] == LibSparseIR.SPIR_COMPUTATION_SUCCESS
@test dlr_u != C_NULL

dlr_uhat_status = Ref{Cint}(-100)
dlr_uhat = LibSparseIR.spir_basis_get_uhat(dlr, dlr_uhat_status)
@test dlr_uhat_status[] == LibSparseIR.SPIR_COMPUTATION_SUCCESS
@test dlr_uhat != C_NULL

# %%
# IR basis functions
ir_u_status = Ref{Cint}(-100)
ir_u = LibSparseIR.spir_basis_get_u(basis, ir_u_status)
@test ir_u_status[] == LibSparseIR.SPIR_COMPUTATION_SUCCESS
@test ir_u != C_NULL

ir_uhat_status = Ref{Cint}(-100)
ir_uhat = LibSparseIR.spir_basis_get_uhat(basis, ir_uhat_status)
@test ir_uhat_status[] == LibSparseIR.SPIR_COMPUTATION_SUCCESS
@test ir_uhat != C_NULL

# %%
function compare_tensors_with_relative_error(a, b, tol)
    diff = abs.(a - b)
    ref = abs.(a)
    maxx_diff = maximum(diff)
    maxx_ref = maximum(ref)
    return maxx_diff <= tol * maxx_ref
end


# %%
function _transform_coeffients(coeffs::AbstractArray{T,N}, basis_eval::AbstractMatrix{U}, target_dim::Integer) where {T,U,N}
    # Move the target dimension to the first position
    coeffs_targetdim0 = movedim(coeffs, 1 + target_dim, 1)
    # Calculate the size of the extra dimensions
    extra_size = prod(coeffs_targetdim0.size[2:end])
    # Create result tensor with correct dimensions
    dims = Vector{Int32}(undef, ndim)
    dims[1] = size(basis_eval, 1)
    dims[2:end] .= size(coeffs_targetdim0)[2:end]
    # Initialize the result
    PromoteType = promote_type(T, U)
    result = Array{PromoteType,N}(undef, dims...)
    # Map tensors to matrices for multiplication
    coeffs_mat = reshape(coeffs_targetdim0, size(coeffs_targetdim0, 1), extra_size)
    result_mat = reshape(result, size(basis_eval, 1), extra_size)
    # Perform the matrix multiplication
    result_mat .= basis_eval * coeffs_mat
    # Move dimension back to original order
    return movedim(result, 1, 1 + target_dim)
end

function _evaluate_basis_functions(::Type{T}, u, x_values) where {T}
    status = Ref{Cint}(-100)
    funcs_size_ref = Ref{Cint}(0)
    status = LibSparseIR.spir_funcs_get_size(u, funcs_size_ref)
    @test status[] == LibSparseIR.SPIR_COMPUTATION_SUCCESS
    funcs_size = funcs_size_ref[]

    u_eval_mat = Matrix{T}(undef, length(x_values), funcs_size)
    for i in eachindex(x_values)
        u_eval = Vector{T}(undef, funcs_size)
        status = LibSparseIR.spir_funcs_eval(u, x_values[i], u_eval)
        @test status[] == LibSparseIR.SPIR_COMPUTATION_SUCCESS
        u_eval_mat[i, :] .= u_eval
    end
    return u_eval_mat
end

function _evaluate_gtau(coeffs::AbstractArray{T,N}, u, target_dim, x_values) where {T,N}
    u_eval_mat = _evaluate_basis_functions(T, u, x_values)
    return _transform_coeffients(coeffs, u_eval_mat, target_dim)
end

function _evaluate_matsubara_basis_functions(uhat, matsubara_indices)
    status = Ref{Cint}(-100)
    funcs_size_ref = Ref{Cint}(0)
    status = LibSparseIR.spir_funcs_get_size(uhat, funcs_size_ref)
    @test status[] == LibSparseIR.SPIR_COMPUTATION_SUCCESS
    funcs_size = funcs_size_ref[]

    uhat_eval_mat = Matrix{ComplexF64}(undef, length(matsubara_indices), funcs_size)
    freq_indices = Int64.(matsubara_indices)
    order = LibSparseIR.SPIR_ORDER_COLUMN_MAJOR
    status = LibSparseIR.spir_funcs_batch_eval_matsu(uhat, order, length(matsubara_indices), freq_indices, uhat_eval_mat)
    @test status[] == LibSparseIR.SPIR_COMPUTATION_SUCCESS
    return uhat_eval_mat
end

function _evaluate_giw(coeffs, uhat, target_dim, matsubara_indices)
    uhat_eval_mat = _evaluate_matsubara_basis_functions(uhat, matsubara_indices)
    result = _transform_coeffients(coeffs, uhat_eval_mat, target_dim)
    return result
end

# %%
num_tau_points_ref = Ref{Cint}(-100)
LibSparseIR.spir_basis_get_n_default_taus(basis, num_tau_points_ref)
num_tau_points = num_tau_points_ref[]
tau_points_org = Vector{T}(undef, num_tau_points)
status = Ref{Cint}(-100)
tau_sampling = LibSparseIR.spir_tau_sampling_new(basis, num_tau_points, tau_points_org, status)
@test status[] == LibSparseIR.SPIR_COMPUTATION_SUCCESS
status = LibSparseIR.spir_sampling_get_npoints(tau_sampling, num_tau_points_ref)
@test status == LibSparseIR.SPIR_COMPUTATION_SUCCESS
num_tau_points = num_tau_points_ref[]
tau_points = Vector{T}(undef, num_tau_points)
status = LibSparseIR.spir_sampling_get_npoints(tau_sampling, tau_points)
@test status == LibSparseIR.SPIR_COMPUTATION_SUCCESS


# %%
# Compare the Greens function at all tau points between IR and DLR
println("Evaluate Greens function at all tau points between IR and DLR")
println("g_IR..")
gtau_from_IR = _evaluate_gtau(g_IR, ir_u, target_dim, tau_points)
println("g_DLR..")
gtau_from_DLR = _evaluate_gtau(coeffs, dlr_u, target_dim, tau_points)
gtau_from_DLR_reconst = _evaluate_gtau(g_DLR_reconst, dlr_u, target_dim, tau_points)

@test compare_tensors_with_relative_error(gtau_from_IR, gtau_from_DLR, tol)
@test compare_tensors_with_relative_error(gtau_from_IR, gtau_from_DLR_reconst, tol)


# %%
num_matsubara_points_org_ref = Ref{Cint}(0)
status = LibSparseIR.spir_basis_get_n_default_matsus(basis, positive_only, num_matsubara_points_org_ref)
@test status == LibSparseIR.SPIR_COMPUTATION_SUCCESS
num_matsubara_points_org = num_matsubara_points_org_ref[]
@test num_matsubara_points_org > 0
matsubara_points_org = Vector{Int64}(undef, num_matsubara_points_org)
status = LibSparseIR.spir_basis_get_default_matsus(basis, positive_only, matsubara_points_org)
@test status == LibSparseIR.SPIR_COMPUTATION_SUCCESS

status = Ref{Cint}(-100)
matsubara_sampling = LibSparseIR.spir_matsu_sampling_new(basis, positive_only, num_matsubara_points_org, matsubara_points_org, status)
@test status[] == LibSparseIR.SPIR_COMPUTATION_SUCCESS

num_matsubara_points_ref = Ref{Cint}(0)
status = LibSparseIR.spir_sampling_get_npoints(matsubara_sampling, num_matsubara_points_ref)
@test status == LibSparseIR.SPIR_COMPUTATION_SUCCESS
num_matsubara_points = num_matsubara_points_ref[]
@test num_matsubara_points > 0
matsubara_points = Vector{Int64}(undef, num_matsubara_points)
status = LibSparseIR.spir_sampling_get_matsus(matsubara_sampling, matsubara_points)
@test status == LibSparseIR.SPIR_COMPUTATION_SUCCESS

# %%
giw_from_IR = _evaluate_giw(g_IR, ir_uhat, target_dim, matsubara_points)
giw_from_DLR = _evaluate_giw(coeffs, dlr_uhat, target_dim, matsubara_points)

@test compare_tensors_with_relative_error(giw_from_IR, giw_from_DLR, tol)

# %%
function _tau_sampling_evaluate(sampling, order, ndim, dims, target_dim, gIR::AbstractArray{<:Real}, gtau::AbstractArray{<:Real})
    status = LibSparseIR.spir_sampling_eval_dd(sampling, order, ndim, dims, target_dim, gIR, gtau)
    @test status == LibSparseIR.SPIR_COMPUTATION_SUCCESS
    return status
end

function _tau_sampling_evaluate(sampling, order, ndim, dims, target_dim, gIR::AbstractArray{<:Complex}, gtau::AbstractArray{<:Complex})
    status = LibSparseIR.spir_sampling_eval_zz(sampling, order, ndim, dims, target_dim, gIR, gtau)
    @test status == LibSparseIR.SPIR_COMPUTATION_SUCCESS
    return status
end

function _tau_sampling_fit(sampling, order, ndim, dims, target_dim, gtau::AbstractArray{<:Real}, gIR::AbstractArray{<:Real})
    status = LibSparseIR.spir_sampling_fit_dd(sampling, order, ndim, dims, target_dim, gtau, gIR)
    @test status == LibSparseIR.SPIR_COMPUTATION_SUCCESS
    return status
end

function _tau_sampling_fit(sampling, order, ndim, dims, target_dim, gtau::AbstractArray{<:Complex}, gIR::AbstractArray{<:Complex})
    status = LibSparseIR.spir_sampling_fit_zz(sampling, order, ndim, dims, target_dim, gtau, gIR)
    @test status == LibSparseIR.SPIR_COMPUTATION_SUCCESS
    return status
end

function _matsubara_sampling_evaluate(sampling, order, ndim, dims, target_dim, gIR::AbstractArray{<:Real}, giw::AbstractArray{<:Complex})
    status = LibSparseIR.spir_sampling_eval_dz(sampling, order, ndim, dims, target_dim, gIR, giw)
    @test status == LibSparseIR.SPIR_COMPUTATION_SUCCESS
    return status
end

function _matsubara_sampling_evaluate(sampling, order, ndim, dims, target_dim, gIR::AbstractArray{<:Complex}, giw::AbstractArray{<:Complex})
    status = LibSparseIR.spir_sampling_eval_zz(sampling, order, ndim, dims, target_dim, gIR, giw)
    @test status == LibSparseIR.SPIR_COMPUTATION_SUCCESS
    return status
end

# %%
dims_matsubara = _get_dims(num_matsubara_points, extra_dims, target_dim, ndim)
dims_IR = _get_dims(basis_size, extra_dims, target_dim, ndim)
dims_tau = _get_dims(num_tau_points, extra_dims, target_dim, ndim)

gIR = Array{T,ndim}(undef, _get_dims(basis_size, extra_dims, target_dim, ndim)...)
gIR2 = Array{T,ndim}(undef, _get_dims(basis_size, extra_dims, target_dim, ndim)...)

gtau = Array{T,ndim}(undef, _get_dims(num_tau_points, extra_dims, target_dim, ndim)...)
giw_reconst = Array{ComplexF64,ndim}(undef, _get_dims(num_matsubara_points, extra_dims, target_dim, ndim)...)

# Matsubara -> IR
begin
    gIR_work = Array{ComplexF64,ndim}(undef, _get_dims(basis_size, extra_dims, target_dim, ndim)...)
    status = LibSparseIR.spir_sampling_fit_zz(
        matsubara_sampling, order, ndim, dims_matsubara, target_dim, giw_from_DLR, gIR_work
    )
    @test status == LibSparseIR.SPIR_COMPUTATION_SUCCESS
    gIR .= real(gIR_work)
end

# IR -> tau
status = _tau_sampling_evaluate(tau_sampling, order, ndim, dims_IR, target_dim, g_IR, gtau)
@test status == LibSparseIR.SPIR_COMPUTATION_SUCCESS

# tau -> IR
status = _tau_sampling_fit(tau_sampling, order, ndim, dims_tau, target_dim, gtau, gIR2)
@test status == LibSparseIR.SPIR_COMPUTATION_SUCCESS

# IR -> Matsubara
status = _matsubara_sampling_evaluate(matsubara_sampling, order, ndim, dims_IR, target_dim, gIR2, giw_reconst)
@test status == LibSparseIR.SPIR_COMPUTATION_SUCCESS

giw_from_IR_reconst = _evaluate_giw(gIR2, ir_uhat, target_dim, matsubara_points)

compare_tensors_with_relative_error(giw_from_DLR, giw_from_IR_reconst, tol)



# %%
