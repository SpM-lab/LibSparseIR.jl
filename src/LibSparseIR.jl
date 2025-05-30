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
       MatsubaraSampling64F, MatsubaraSampling64B, TauSampling64, sampling_points,
       basis

include("freq.jl")
include("abstract.jl")
include("kernel.jl")
include("sve.jl")

end # module LibSparseIR
