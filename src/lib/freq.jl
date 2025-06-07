"""
    Statistics(zeta)

Abstract type for quantum statistics (fermionic/bosonic/etc.)
"""
abstract type Statistics end

function Statistics(zeta::Integer)
    if isone(zeta)
        return Fermionic()
    elseif iszero(zeta)
        return Bosonic()
    else
        throw(DomainError(zeta, "does not correspond to known statistics"))
    end
end


"""
Fermionic statistics.
"""
struct Fermionic <: Statistics end

"""
Bosonic statistics.
"""
struct Bosonic <: Statistics end


# Convert Julia statistics to C API constants
_statistics_to_c(::Type{Fermionic}) = SPIR_STATISTICS_FERMIONIC
_statistics_to_c(::Type{Bosonic}) = SPIR_STATISTICS_BOSONIC
_statistics_from_c(s::Cint) = s == SPIR_STATISTICS_FERMIONIC ? Fermionic() : Bosonic()


"""
    MatsubaraFreq(n)

Prefactor `n` of the Matsubara frequency `ω = n*π/β`

Struct representing the Matsubara frequency ω entering the Fourier transform of
a propagator G(τ) on imaginary time τ to its Matsubara equivalent Ĝ(iω) on the
imaginary-frequency axis:

            β
    Ĝ(iω) = ∫  dτ exp(iωτ) G(τ)      with    ω = n π/β,
            0

where β is inverse temperature and by convention we include the imaginary unit
in the frequency argument, i.e, Ĝ(iω). The frequencies depend on the
statistics of the propagator, i.e., we have that:

    G(τ - β) = ± G(τ)

where + is for bosons and - is for fermions. The frequencies are restricted
accordingly.

  - Bosonic frequency (`S == Fermionic`): `n` even (periodic in β)
  - Fermionic frequency (`S == Bosonic`): `n` odd (anti-periodic in β)
"""
struct MatsubaraFreq{S<:Statistics} <: Number
    n::Int

    MatsubaraFreq(stat::Statistics, n::Integer) = new{typeof(stat)}(n)

    function MatsubaraFreq{S}(n::Integer) where {S<:Statistics}
        allowed(S, n) || throw(DomainError(n, "Frequency $(n)π/β is not $S"))
        return new{S}(n)
    end
end

const BosonicFreq   = MatsubaraFreq{Bosonic}
const FermionicFreq = MatsubaraFreq{Fermionic}

MatsubaraFreq(n::Integer) = MatsubaraFreq(Statistics(mod(n, 2)), n)
