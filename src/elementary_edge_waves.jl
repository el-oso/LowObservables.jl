"""
3-D electromagnetic elementary edge-wave (EEW) directivity coefficients F⁽¹⁾ and G⁽¹⁾.

These are the *spherical* (ϑ,φ) components of the EEW vectors, valid for **all** observation
directions — not just on the Keller diffraction cone.  On the cone (ϑ = π − γ₀) they reduce
to the 2-D fringe functions f⁽¹⁾/g⁽¹⁾ via Eqs. 7.148 and 7.151.  G_ϑ⁽¹⁾ (G_theta) is an
exclusively electromagnetic polarisation-coupling term with no acoustic analog (Eq. 7.153).

Source: Ufimtsev (2014), §7.8.3 and *Appendix to Section 7.8.3*,
MATLAB listing `FG.m`, pp. 443–445.  Equations 7.75–7.94, 7.106–7.109, 7.139–7.141, 7.148–7.153.

# ponytail: TypeContracts @contract skipped — standalone function, same rationale as fringe.jl.
# ponytail: off-cone → on-cone limit test deferred; on-cone branch is validated by construction
# (it calls fringe_directivity directly). Add a limit test if an off-cone oracle emerges.
"""

"""
    eew_directivity(φ, φ₀, γ₀, ϑ, α) -> (; F_theta, F_phi, G_theta, G_phi)

Spherical (ϑ,φ) components of the electromagnetic EEW directivity coefficients F⁽¹⁾ and G⁽¹⁾
for a perfectly-conducting wedge.  Faithful port of `FG.m`, Ufimtsev (2014), pp. 443–445.

# Arguments
- `φ`  : observation azimuth [rad], `0 ≤ φ ≤ α`
- `φ₀` : incidence azimuth [rad], `0 < φ₀ < π`
- `γ₀` : diffraction-cone half-angle [rad] (incidence polar angle from the edge)
- `ϑ`  : observation polar angle [rad], `0 < ϑ < π`
- `α`  : exterior wedge angle [rad]; `α > π` for a convex edge (e.g. `3π/2` for 270°)

# Returns NamedTuple
- `F_theta` : ϑ-component of F⁽¹⁾ (soft / E-pol)
- `F_phi`   : φ-component of F⁽¹⁾ (always 0 by symmetry — Eq. 7.139)
- `G_theta` : ϑ-component of G⁽¹⁾ (hard / H-pol; includes EM polarisation-coupling)
- `G_phi`   : φ-component of G⁽¹⁾ (hard / H-pol)

## On the Keller cone (`ϑ = π − γ₀`):
  `F_theta = −f1/sin(γ₀)`,  `G_phi = g1/sin(γ₀)`  (Eqs. 7.148, 7.151)
  `G_theta = (ε(φ₀) − ε(α−φ₀))·cot(γ₀)`  (Eq. 7.153, no acoustic analog)
  where ε(x) = 1 for x ∈ (0,π), else 0.

## Note on singularities
F⁽¹⁾ and G⁽¹⁾ are singular at the two grazing illumination directions
`φ₀ = π` (with `φ = 0`) and `φ₀ = α − π` (with `φ = α`).
Shadow/reflection boundaries off the cone are handled by implicit cancellation.

Generic over float type: `Float32` inputs → `Float32` outputs.
"""
function eew_directivity(φ::Real, φ₀::Real, γ₀::Real, ϑ::Real, α::Real)
    T = promote_type(typeof(float(φ)), typeof(float(φ₀)), typeof(float(γ₀)),
                     typeof(float(ϑ)), typeof(float(α)))
    _eew_directivity(T(φ), T(φ₀), T(γ₀), T(ϑ), T(α))
end

function _eew_directivity(φ::T, φ₀::T, γ₀::T, ϑ::T, α::T) where {T<:AbstractFloat}
    # ── On-cone branch: Vtheta == pi - gamma0 ────────────────────────────────
    # Reuse Phase-5a fringe_directivity; on-cone formulas are Eqs. 7.148–7.153.
    if abs(ϑ - (T(π) - γ₀)) < T(1e-10)
        fg = fringe_directivity(φ, φ₀, α)
        return (;
            F_theta = -fg.f1 / sin(γ₀),                               # Eq.(7.148)
            F_phi   = zero(T),                                          # Eq.(7.150)
            G_theta = (_eps_x(φ₀) - _eps_x(α - φ₀)) * cot(γ₀),      # Eq.(7.153)
            G_phi   = fg.g1 / sin(γ₀),                                 # Eq.(7.151)
        )
    end

    # ── Off-cone general case (Eqs. 7.75–7.141) ──────────────────────────────
    # β₁, β₂: polar angles from the edge to the observation point, projected on each face.
    β₁ = acos(clamp(sin(γ₀)*sin(ϑ)*cos(φ)       - cos(γ₀)*cos(ϑ), T(-1), T(1)))  # Eq.(7.75)
    β₂ = acos(clamp(sin(γ₀)*sin(ϑ)*cos(α - φ)   - cos(γ₀)*cos(ϑ), T(-1), T(1)))  # Eq.(7.95)

    # σ₁, σ₂: effective azimuths on the diffraction cone (auxiliary angles).
    # sigma12(β, γ₀) = arccos((cosβ − cos²γ₀)/sin²γ₀)
    σ₁ = _sigma12(β₁, γ₀)
    σ₂ = _sigma12(β₂, γ₀)

    s2 = sin(γ₀)^2
    φ₂ = α - φ₀   # complement incidence azimuth (face 2 reference)

    # Face-1 building blocks: Eqs.(7.85)–(7.88)
    Ut1 = (T(π) / (2α*s2)) * (cot(T(π)*(σ₁ + φ₀)/(2α)) - cot(T(π)*(σ₁ - φ₀)/(2α)))  # Eq.(7.85)
    U01 = -_eps_x(φ₀) * sin(φ₀) / (cos(β₁) - cos(γ₀)^2 + s2*cos(φ₀))                  # Eq.(7.86)
    Vt1 = (T(π) / (2α*s2*sin(σ₁))) * (cot(T(π)*(σ₁ + φ₀)/(2α)) + cot(T(π)*(σ₁ - φ₀)/(2α)))  # Eq.(7.87)
    V01 = _eps_x(φ₀) / (cos(β₁) - cos(γ₀)^2 + s2*cos(φ₀))                              # Eq.(7.88)

    # Face-2 building blocks: symmetric with φ₂ = α − φ₀ (analogous to Eqs. 7.85–7.88)
    Ut2 = (T(π) / (2α*s2)) * (cot(T(π)*(σ₂ + φ₂)/(2α)) - cot(T(π)*(σ₂ - φ₂)/(2α)))
    U02 = -_eps_x(φ₂) * sin(φ₂) / (cos(β₂) - cos(γ₀)^2 + s2*cos(φ₂))
    Vt2 = (T(π) / (2α*s2*sin(σ₂))) * (cot(T(π)*(σ₂ + φ₂)/(2α)) + cot(T(π)*(σ₂ - φ₂)/(2α)))
    V02 = _eps_x(φ₂) / (cos(β₂) - cos(γ₀)^2 + s2*cos(φ₂))

    # Fringe = total − PO: Eqs.(7.83),(7.84),(7.93),(7.94)
    U1 = Ut1 - U01;  V1 = Vt1 - V01
    U2 = Ut2 - U02;  V2 = Vt2 - V02

    # Shadow-boundary regularisation (removable pole when σ = φ₀):
    # Eqs.(7.106),(7.108),(7.109)
    if abs(σ₁ - φ₀) < T(1e-10)
        U1 = T(π)*cot(T(π)*φ₀/α) / (2α*s2) - cot(φ₀) / (2s2)   # Eq.(7.106)
        V1 = U1 / sin(φ₀)                                          # Eq.(7.109)
    end
    if abs(σ₂ - φ₂) < T(1e-10)
        U2 = T(π)*cot(T(π)*φ₂/α) / (2α*s2) - cot(φ₂) / (2s2)   # Eq.(7.108)
        V2 = U2 / sin(φ₂)                                          # Eq.(7.109)
    end

    # Spherical components: Eqs.(7.139)–(7.141)
    F_theta = (U1 + U2) * sin(ϑ)                                   # Eq.(7.139)
    F_phi   = zero(T)                                               # Eq.(7.139)
    G_theta = sin(ϑ)*cos(γ₀)*(_eps_x(φ₀) - _eps_x(φ₂)) / s2 +
              (sin(γ₀)*cos(ϑ)*cos(φ)       - cos(γ₀)*sin(ϑ)*cos(σ₁)) * V1 -
              (sin(γ₀)*cos(ϑ)*cos(α - φ)   - cos(γ₀)*sin(ϑ)*cos(σ₂)) * V2  # Eq.(7.140)
    G_phi   = -(V1*sin(φ) + V2*sin(α - φ)) * sin(γ₀)              # Eq.(7.141)

    return (; F_theta, F_phi, G_theta, G_phi)
end

# ε(x) = 1 for x ∈ (0, π), else 0.  Illuminated-face indicator (Heaviside-step pair).
_eps_x(x::T) where {T<:AbstractFloat} = T(0) < x < T(π) ? one(T) : zero(T)
_eps_x(x::Real) = _eps_x(float(x))

# Effective azimuth on the diffraction cone σ(β,γ) — faithful port of sigma12.m,
# Ufimtsev (2014) p.441, Eqs.(7.76)/(7.77).  REAL on the near-cone branch (β ≤ βₖ),
# COMPLEX on the wide-angle branch (β > βₖ).  Always returned as Complex{T}; the
# imaginary part is 0 on the near-cone branch (continuous reduction to the 2-D σ = π−φ).
# ponytail: the real↔complex branches do NOT meet at β = βₖ — the published sigma12 is
# discontinuous on that fold (β₂ = 2γ₀, reached e.g. when φ = α/2).  This is an artifact of
# the book's algorithm, faithfully reproduced; smoothing it would deviate from the source.
function _sigma12(β::T, γ₀::T) where {T<:AbstractFloat}
    βk = (T(0) <= γ₀ <= T(π)/2) ? 2γ₀ : 2*(T(π) - γ₀)
    if β <= βk
        return complex(T(π) - acos(clamp((cos(β) - cos(γ₀)^2) / sin(γ₀)^2, T(-1), T(1))))  # Eq.(7.76)
    else
        return im*log(complex(cos(γ₀)^2 - cos(β) + sqrt(complex(cos(β)^2 - sin(γ₀)^4)))) -
               2im*log(sin(γ₀))                                                              # Eq.(7.77)
    end
end
