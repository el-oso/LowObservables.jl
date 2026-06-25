"""
Mie-series RCS for a perfectly-conducting (PEC) sphere.

References:
- Bohren & Huffman, *Absorption and Scattering of Light by Small Particles* (1983), §4.4
- Wiscombe (1980), "Improved Mie Scattering Algorithms", Appl. Opt. 19, 1505
"""

using SpecialFunctions: besselj, bessely

# ── Riccati-Bessel functions ──────────────────────────────────────────────────

"""
    riccati_j(n, x) -> (ψ_n, ψ_{n-1})

Riccati-Bessel function ψ_n(x) = x·j_n(x) and ψ_{n-1}(x) = x·j_{n-1}(x),
computed from `besselj` at half-integer orders.
"""
function riccati_j(n::Int, x::T) where {T<:Real}
    c = sqrt(T(π) * x / 2)
    ψn  = c * besselj(n + one(T)/2, x)
    ψnm1 = c * besselj(n - one(T)/2, x)
    return ψn, ψnm1
end

"""
    riccati_h(n, x) -> (ξ_n, ξ_{n-1})

Riccati-Hankel function ξ_n(x) = x·h_n^(1)(x) and ξ_{n-1}(x), both complex.
h_n^(1) = j_n + i·y_n, so ξ_n = x·j_n + i·x·y_n.
"""
function riccati_h(n::Int, x::T) where {T<:Real}
    c = sqrt(T(π) * x / 2)
    jn  = c * besselj(n + one(T)/2, x)
    yn  = c * bessely(n + one(T)/2, x)
    jnm1 = c * besselj(n - one(T)/2, x)
    ynm1 = c * bessely(n - one(T)/2, x)
    ξn   = jn   + im * yn
    ξnm1 = jnm1 + im * ynm1
    return ξn, ξnm1
end

# ── Mie coefficients for a PEC sphere ────────────────────────────────────────

"""
    pec_mie_an_bn(n, x) -> (a_n, b_n)

Mie scattering coefficients for a PEC sphere at size parameter x = ka.

For PEC (E·n̂ = 0 at surface):
  a_n = ψ_n(x) / ξ_n(x)       (TE / electric Mie coefficient)
  b_n = ψ_n'(x) / ξ_n'(x)     (TM / magnetic Mie coefficient)

Derivatives use the Riccati-Bessel recurrence:
  f_n'(x) = f_{n-1}(x) - (n/x)·f_n(x)
"""
function pec_mie_an_bn(n::Int, x::T) where {T<:Real}
    ψn, ψnm1 = riccati_j(n, x)
    ξn, ξnm1 = riccati_h(n, x)

    # derivatives via recurrence
    nx = T(n) / x
    dψn = ψnm1 - nx * ψn
    dξn = ξnm1 - nx * ξn

    an = ψn / ξn
    bn = dψn / dξn
    return an, bn
end

# ── Wiscombe truncation criterion ─────────────────────────────────────────────

"""
    wiscombe_nmax(x) -> Int

Number of Mie terms to sum for size parameter x = ka.
Wiscombe (1980): N = x + 4x^(1/3) + 2, minimum 5.
"""
wiscombe_nmax(x::Real) = max(5, ceil(Int, x + 4 * x^(one(x)/3) + 2))

# ── Public API ────────────────────────────────────────────────────────────────

"""
    mie_rcs_sphere(a, λ) -> σ

Monostatic (backscatter) radar cross-section [m²] of a perfectly-conducting
sphere of radius `a` [m] at wavelength `λ` [m], via the exact Mie series.

The formula (Bohren & Huffman §4.4):
  σ = (λ²/π)|S(π)|²
  S(π) = (1/2) Σ_n (-1)^n (2n+1)(b_n - a_n)

where a_n, b_n are the PEC Mie coefficients and the sum runs to N≈ka+4(ka)^(1/3)+2.

Generic over the float type of `a` and `λ`: pass `Float32` inputs → `Float32` RCS.
"""
function mie_rcs_sphere(a::Real, λ::Real)
    T = promote_type(typeof(float(a)), typeof(float(λ)))
    af = T(a)
    λf = T(λ)
    k  = 2 * T(π) / λf
    x  = k * af

    N = wiscombe_nmax(x)

    S = zero(Complex{T})
    for n in 1:N
        an, bn = pec_mie_an_bn(n, x)
        sign_n = iseven(n) ? one(T) : -one(T)     # (-1)^n
        S += sign_n * T(2n + 1) * (bn - an)
    end

    # backscatter RCS: σ = (λ²/π)|S(π)|²  with S(π) = S_sum/2
    σ = (λf^2 / T(π)) * abs2(S / 2)
    return σ
end

"""
    rcs_vs_size(a, λs) -> Vector

Monostatic RCS [m²] for a PEC sphere of radius `a` over a vector of wavelengths `λs`.
Convenience wrapper around `mie_rcs_sphere`.
"""
rcs_vs_size(a::Real, λs::AbstractVector{<:Real}) = [mie_rcs_sphere(a, λ) for λ in λs]

"""
    optical_limit_rcs(a) -> σ_opt

Geometric-optics reference cross-section πa² [m²] for a sphere of radius `a`.
"""
optical_limit_rcs(a::Real) = oftype(float(a), π) * float(a)^2

# ── MieSphere model type ──────────────────────────────────────────────────────

"""
    MieSphere <: AbstractRCSModel

PEC sphere of radius `radius` [m], using the exact Mie series.
Radius is stored as Float64; for other float types call `mie_rcs_sphere` directly.

```julia
m = MieSphere(0.1)          # 10 cm sphere
rcs(m, 0.03)                # RCS at λ = 3 cm (≈ 10 GHz)
```

# ponytail: non-parametric so TypeContracts @verify can check the return type
# statically (Base.return_types needs a concrete type, not a UnionAll).
# Users needing Float32/AD use mie_rcs_sphere(Float32(a), Float32(λ)) directly.
"""
struct MieSphere <: AbstractRCSModel
    radius::Float64
end

"""
    rcs(m::MieSphere, λ) -> Float64

Monostatic Mie RCS [m²] for the sphere at wavelength λ [m].
"""
rcs(m::MieSphere, λ::Real)::Float64 = mie_rcs_sphere(m.radius, Float64(λ))

"""
    rcs_vs_wavelength(m::MieSphere, λs) -> Vector

Mie RCS [m²] over a vector of wavelengths.
"""
rcs_vs_wavelength(m::MieSphere, λs::AbstractVector{<:Real}) = rcs_vs_size(m.radius, λs)

@verify MieSphere
