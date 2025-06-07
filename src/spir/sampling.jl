# Convenience property accessors (similar to SparseIR.jl)
Base.getproperty(s::TauSampling, p::Symbol) = p === :tau ? sampling_points(s) : getfield(s, p)
Base.getproperty(s::MatsubaraSampling, p::Symbol) = p === :wn ? sampling_points(s) : getfield(s, p)
