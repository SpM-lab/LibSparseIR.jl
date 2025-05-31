module LibSparseIR

include("C_API.jl") # LibSparseIR_C_API
using .C_API

import LinearAlgebra

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

include("freq.jl")
include("abstract.jl")
include("kernel.jl")
include("sve.jl")
include("basis.jl")
include("sampling.jl")
include("dlr.jl")

end # module LibSparseIR
