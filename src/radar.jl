"""
Radar Range Equation — Phase 6.

Monostatic single-antenna radar equation, detection range, radar horizon,
and detection-range sweep over aspect angle.

# Physics

SNR at slant range R:
  SNR = (Pt · G² · λ² · σ) / ((4π)³ · R⁴ · kB · T0 · B · F · L)

Detection range (max R for SNR_min):
  R_max = ((Pt · G² · λ² · σ) / ((4π)³ · kB · T0 · B · F · L · SNR_min))^(1/4)

KEY INSIGHT: R_max ∝ σ^(1/4) — a ×10 RCS reduction cuts R_max only ×10^0.25 ≈ 1.78.

Radar horizon (standard 4/3-earth refraction):
  R_hor ≈ 4.12 · (√h_target + √h_radar)  [km, h in metres]

# Reference
Knott, Schaeffer, Tuley, *Radar Cross Section*, 2nd ed. (2004) §2.5, §14.2

# ponytail: TypeContracts @contract skipped — no abstract-type interface;
#   see bodies_of_revolution.jl for the convention.
# ponytail: Swerling fluctuation models, atmospheric attenuation, pulse integration
#   deferred; add when the use-case warrants it.
"""

# ── Physical constants ─────────────────────────────────────────────────────────

const _kB = 1.380649e-23   # Boltzmann constant [J/K]
const _T0  = 290.0          # standard noise temperature [K]
const _c   = 3e8            # speed of light [m/s]

# ── dB helpers ────────────────────────────────────────────────────────────────

"""
    db2lin(x) → linear power ratio

Convert from dB to linear (power). Example: `db2lin(10) == 10.0`.
"""
db2lin(x::Real) = 10^(x / 10)

"""
    lin2db(x) → dB

Convert from linear power ratio to dB. Example: `lin2db(100) == 20.0`.
"""
lin2db(x::Real) = 10 * log10(x)

# ── RadarSystem ───────────────────────────────────────────────────────────────

"""
    RadarSystem{T<:AbstractFloat}

Parameters for a monostatic pulsed-Doppler radar.

# Fields
- `Pt`      — peak transmit power [W]
- `G`       — antenna gain (linear; `db2lin(G_dB)`)
- `freq`    — centre frequency [Hz]  (λ = c/freq)
- `B`       — noise bandwidth [Hz]
- `F`       — noise figure (linear; `db2lin(NF_dB)`)
- `L`       — system losses (linear; `db2lin(L_dB)`)
- `SNR_min` — minimum detection SNR (linear)

# Constructor

All fields have defaults modelling an X-band air-defence radar (Pt=1 MW, G=40 dB,
f=10 GHz, B=1 MHz, NF≈4.8 dB, L≈6 dB, SNR_min≈11 dB):

    RadarSystem(; Pt=1e6, G=1e4, freq=10e9, B=1e6, F=3.0, L=4.0, SNR_min=13.0)

Pass `Float32` values for a `RadarSystem{Float32}`.
"""
struct RadarSystem{T<:AbstractFloat}
    Pt      :: T
    G       :: T
    freq    :: T
    B       :: T
    F       :: T
    L       :: T
    SNR_min :: T
end

function RadarSystem(;
    Pt      :: Real = 1e6,
    G       :: Real = 1e4,
    freq    :: Real = 10e9,
    B       :: Real = 1e6,
    F       :: Real = 3.0,
    L       :: Real = 4.0,
    SNR_min :: Real = 13.0,
)
    T = promote_type(
        typeof(float(Pt)), typeof(float(G)), typeof(float(freq)),
        typeof(float(B)),  typeof(float(F)), typeof(float(L)),
        typeof(float(SNR_min)),
    )
    RadarSystem{T}(T(Pt), T(G), T(freq), T(B), T(F), T(L), T(SNR_min))
end

# ── Core equations ─────────────────────────────────────────────────────────────

"""
    radar_snr(radar, σ, R) → SNR (linear)

Monostatic radar SNR at slant range `R` [m] for a target with RCS `σ` [m²].

    SNR = (Pt·G²·λ²·σ) / ((4π)³·R⁴·kB·T0·B·F·L)

Convert to dB with `lin2db(radar_snr(...))`.  SNR ∝ 1/R⁴: doubling range
drops SNR by 16×.
"""
function radar_snr(radar::RadarSystem{T}, σ::Real, R::Real) where {T}
    λ   = T(_c) / radar.freq
    num = radar.Pt * radar.G^2 * λ^2 * T(σ)
    den = (4 * T(π))^3 * T(R)^4 * T(_kB) * T(_T0) * radar.B * radar.F * radar.L
    return num / den
end

"""
    detection_range(radar, σ) → R_max [m]

Maximum detection range for a target with RCS `σ` [m²], at the radar's minimum
detectable SNR.

    R_max = ((Pt·G²·λ²·σ) / ((4π)³·kB·T0·B·F·L·SNR_min))^(1/4)

The **R⁴ dependence** (equivalently σ^(1/4) scaling) is the central result:
  - ×10 in RCS  → ×10^0.25 ≈ 1.78 in detection range (not ×10)
  - ×10000 in RCS → ×10 in range
This is why stealth technology requires enormous RCS reductions for modest
range benefits, and why shaping + RAM must be extremely effective.
"""
function detection_range(radar::RadarSystem{T}, σ::Real) where {T}
    λ   = T(_c) / radar.freq
    num = radar.Pt * radar.G^2 * λ^2 * T(σ)
    den = (4 * T(π))^3 * T(_kB) * T(_T0) * radar.B * radar.F * radar.L * radar.SNR_min
    return (num / den)^T(0.25)
end

"""
    radar_horizon(h_target; h_radar=0) → R_horizon [m]

Line-of-sight radar horizon under standard 4/3-earth refraction.

    R_hor ≈ 4.12·(√h_target + √h_radar)  [km]  (h in metres)

**Returns metres.**  A target at slant range R > R_horizon is undetectable
regardless of RCS or radar power.  For a ground-based radar (`h_radar=0`) and
a target at 10 km altitude: R_hor ≈ 412 km.
"""
function radar_horizon(h_target::Real; h_radar::Real = zero(h_target))
    T = promote_type(typeof(float(h_target)), typeof(float(h_radar)))
    return T(4120) * (sqrt(T(h_target)) + sqrt(T(h_radar)))   # [m]
end

"""
    detectable(radar, σ, R; h_target, h_radar=0) → Bool

True iff the target satisfies both the range-equation gate and the LOS horizon gate:

    R ≤ detection_range(radar, σ)  AND  R ≤ radar_horizon(h_target; h_radar)

A target beyond either limit is undetectable.
"""
function detectable(
    radar    :: RadarSystem,
    σ        :: Real,
    R        :: Real;
    h_target :: Real,
    h_radar  :: Real = 0.0,
)
    return R ≤ detection_range(radar, σ) && R ≤ radar_horizon(h_target; h_radar)
end

# ── Aspect sweep ──────────────────────────────────────────────────────────────

"""
    detection_range_vs_aspect(mesh, dirs, ei; radar, k) → Vector{T}

Maximum detection range [m] for each incidence direction in `dirs` (3×N matrix of
unit vectors).  RCS at each aspect is computed by `ptd_rcs_monostatic` (PO + PTD
edge corrections) and fed into `detection_range`.

Generic over the float type of `mesh` (Float32 and Float64 supported).

# Arguments
- `mesh`  — PEC surface mesh
- `dirs`  — 3×N matrix of unit incident direction vectors k̂_i (each column)
- `ei`    — incident E-field polarisation (unit vector, ⊥ ki for all columns)
- `radar` — `RadarSystem`
- `k`     — wavenumber 2π/λ [rad/m]
"""
function detection_range_vs_aspect(
    mesh  :: TriMesh{T},
    dirs  :: AbstractMatrix{<:Real},
    ei    :: AbstractVector{<:Real};
    radar :: RadarSystem,
    k     :: Real,
) where {T<:AbstractFloat}
    n      = size(dirs, 2)
    result = Vector{T}(undef, n)
    for j in 1:n
        ki_j = (T(dirs[1,j]), T(dirs[2,j]), T(dirs[3,j]))
        σ    = ptd_rcs_monostatic(mesh, ki_j, ei; k)
        result[j] = T(detection_range(radar, σ))
    end
    return result
end
