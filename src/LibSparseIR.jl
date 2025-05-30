module LibSparseIR

using LinearAlgebra: LinearAlgebra

include("C_API.jl") # LibSparseIR_C_API
using .C_API

end # module LibSparseIR
