module LibSparseIR

include("C_API.jl") # LibSparseIR_C_API
using .C_API

import LinearAlgebra
using QuadGK: quadgk

export Fermionic, Bosonic
export MatsubaraFreq, BosonicFreq, FermionicFreq, pioverbeta
export FiniteTempBasis, FiniteTempBasisSet
export DiscreteLehmannRepresentation
export overlap
export LogisticKernel, RegularizedBoseKernel
export AugmentedBasis, TauConst, TauLinear, MatsubaraConst
export TauSampling, MatsubaraSampling, evaluate, fit, evaluate!, fit!,
       sampling_points, basis, npoints
export from_IR, to_IR, npoles, get_poles, default_omega_sampling_points

include("lib/freq.jl")
include("lib/abstract.jl")
include("lib/kernel.jl")
include("lib/sve.jl")
include("lib/poly.jl")
include("lib/basis.jl")
include("lib/sampling.jl")
include("lib/dlr.jl")

include("spir/abstract.jl")
include("spir/freq.jl")
include("spir/kernel.jl")
include("spir/poly.jl")
include("spir/basis.jl")
include("spir/sampling.jl")
include("spir/dlr.jl")
include("spir/basis_set.jl")
include("spir/augment.jl")

end # module LibSparseIR
