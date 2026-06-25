"""
Closed-form monostatic RCS for canonical bodies of revolution (axial backscattering).

These are the *validation answer-key* formulas — exact (within the PTD approximation)
results that the Physical Optics solver will be validated against in later phases.

Source: Ufimtsev, *Fundamentals of the Physical Theory of Diffraction*, 2nd ed. (2014),
Appendix to Chapter 6, pp. 431–437.  MATLAB files: SCSn_general.m, Fig_610.m – Fig_619.m.

## Normalization conventions

`scsn_general` and its wrappers return the **normalized SCS = σ/(πa²)** where `a` is
the rim radius.  In the Physical-Optics limit:
  - Flat disk at normal incidence: PO = (ka)²
  - Sphere front hemisphere (w = 0): PO = 1.0 exactly (σ_PO = πa² for all ka)

`cone_rcs` returns **σ·k²/π** (no fixed reference area; the cone has no natural radius):
  - abs(Ish0)² is the "normalized SCS" as plotted in Figs. 6.10–6.12.

Both functions are generic over the real float type: pass `Float32` inputs → `Float32` output.

# ponytail: TypeContracts @contract/@verify applies to abstract-type interface contracts
# (AbstractRCSModel).  Standalone functions scsn_general/cone_rcs are not interface
# implementations; their contract is expressed through docstrings + Julia return types.
# Add a CanonicalBOR abstract type + @contract if multiple BOR implementations emerge.
"""

# ── scsn_general ──────────────────────────────────────────────────────────────

"""
    scsn_general(ka, kl, w, Omega, shape) -> (; PO, PTD_soft, PTD_hard)

Normalized monostatic backscatter cross-section (axial incidence) for a
paraboloid or spherical body of revolution.

Port of `SCSn_general.m`, Ufimtsev (2014) Appendix to Chapter 6, p. 437.

# Arguments
- `ka`    : k·a  — wavenumber × rim radius (dimensionless ≥ 0)
- `kl`    : k·l  — wavenumber × axial body length (dimensionless ≥ 0)
- `w`     : rim (half-angle) [rad]; `w = 0` → sphere front hemisphere, `w = π/2` → flat disk
- `Omega` : exterior (shadow-boundary) half-angle [rad]; typically `π/2`
- `shape` : `:Paraboloid` or `:Spherical`

# Returns NamedTuple
- `PO`       : σ/(πa²), Physical-optics approximation (same for soft & hard BC)
- `PTD_soft` : σ/(πa²), PTD (Dirichlet / soft BC)
- `PTD_hard` : σ/(πa²), PTD (Neumann / hard BC)

The special case `w = π/2` (flat-disk limit) is triggered by passing exactly
the floating-point value `T(π)/2` where T is the promoted float type.
The `kl` argument is ignored in that branch (flat-base formula is kl-independent).

Generic over float type: pass `Float32` inputs → `Float32` output.
"""
function scsn_general(ka::Real, kl::Real, w::Real, Omega::Real, shape::Symbol)
    shape ∈ (:Paraboloid, :Spherical) ||
        throw(ArgumentError("shape must be :Paraboloid or :Spherical, got :$shape"))

    T   = promote_type(typeof(float(ka)), typeof(float(kl)),
                       typeof(float(w)),  typeof(float(Omega)))
    kaf = T(ka); klf = T(kl); wf = T(w); Omf = T(Omega)

    n   = 1 + (wf + Omf) / T(π)     # exterior-angle index n (Ufimtsev §6, App. p. 437)

    # Complex scattering amplitudes Ish0 (PO), Is (PTD soft), Ih (PTD hard).
    # kR is computed only in the else branch to avoid 1/cos(π/2) = Inf for the Spherical shape.
    Ish0, Is, Ih = if wf == T(π) / 2
        # ── Flat-base limit  (SCSn_general.m, w == pi/2 branch, p. 437) ──────
        # Body degenerates to a flat circular disk (paraboloid as kl→0, or spherical
        # segment at 90° rim angle).
        sn  = sin(T(π) / n)
        cn  = cos(T(π) / n)
        ctn = cn / sn                              # cot(π/n)
        i0  = Complex{T}(kaf * kaf)               # Ish0 = ka²  (PO aperture integral)
        is  = kaf * (im * kaf + (1/n) * ctn + (2/n) * sn / (cn - 1))
        ih  = kaf * (im * kaf + (1/n) * ctn - (2/n) * sn / (cn - 1))
        i0, is, ih
    else
        # ── General case  (SCSn_general.m, else branch, p. 437) ──────────────
        # kR = k·R, where R is the tip curvature radius at ρ = 0 (axis).
        #   Paraboloid: R = a·tan(w)   (focal-point curvature of the parabolic apex)
        #   Spherical:  R = a/cos(w)   (sphere radius; larger than rim radius a for w > 0)
        kR = shape === :Paraboloid ? kaf * tan(wf) : kaf / cos(wf)
        e  = exp(im * 2 * klf)                    # phase factor e^{i·2kl}
        sn = sin(T(π) / n)
        cn = cos(T(π) / n)
        d1 = cn - one(T)                           # cos(π/n) − 1
        d2 = cn - cos(2 * wf / n)                 # cos(π/n) − cos(2w/n)
        i0 = kR - kaf * tan(wf) * e
        is = kR - kaf * (2/n) * sn * (1/d1 - 1/d2) * e
        ih = kR + kaf * (2/n) * sn * (1/d1 + 1/d2) * e
        i0, is, ih
    end

    ka2 = kaf^2
    return (;
        PO       = abs2(Ish0) / ka2,   # σ/(πa²), PO approximation
        PTD_soft = abs2(Is)   / ka2,   # σ/(πa²), PTD soft (Dirichlet) BC
        PTD_hard = abs2(Ih)   / ka2,   # σ/(πa²), PTD hard (Neumann) BC
    )
end

# ── cone_rcs ──────────────────────────────────────────────────────────────────

"""
    cone_rcs(ka, kl, w, Omega) -> (; PO, PTD_soft, PTD_hard)

Normalized monostatic backscatter cross-section (axial incidence) for a finite cone.

Port of the cone driver code in `Fig_610.m`–`Fig_612.m`,
Ufimtsev (2014) Appendix to Chapter 6, pp. 432–434.

# Arguments
- `ka`    : k·a  — wavenumber × base radius
- `kl`    : k·l  — wavenumber × axial (slant) length; for a right cone `ka = kl·sin(w)`
- `w`     : cone half-angle from the axis [rad]
- `Omega` : exterior (shadow-boundary) half-angle [rad]; `π/2` for a solid cone

# Returns NamedTuple
- `PO`       : σ·k²/π, Physical-optics approximation
- `PTD_soft` : σ·k²/π, PTD soft (Dirichlet) BC
- `PTD_hard` : σ·k²/π, PTD hard (Neumann) BC

The normalization σ·k²/π differs from `scsn_general`'s σ/(πa²) because a cone has
no intrinsic reference area.

Note: the original MATLAB code has an unreachable `elseif (w==pi/2 && omeg==pi/2)` branch.
This implementation corrects the condition order: the more specific double-special-case
is checked first, making both branches reachable.
"""
function cone_rcs(ka::Real, kl::Real, w::Real, Omega::Real)
    T   = promote_type(typeof(float(ka)), typeof(float(kl)),
                       typeof(float(w)),  typeof(float(Omega)))
    kaf = T(ka); klf = T(kl); wf = T(w); Omf = T(Omega)

    n = 1 + (wf + Omf) / T(π)       # exterior-angle index

    # Condition order: more specific case (w=π/2 AND Ω=π/2) checked first —
    # this corrects the unreachable-branch bug present in the original MATLAB.
    Ish0, Is, Ih = if wf == T(π) / 2 && Omf == T(π) / 2
        # ── Knife-edge flat disk: w = π/2, Ω = π/2  (Fig_610.m elseif, p. 432) ──
        i0 = Complex{T}(kaf)
        is = im * kaf - one(T)
        ih = im * kaf + one(T)
        i0, is, ih
    elseif wf == T(π) / 2
        # ── Flat-base cone: w = π/2, Ω ≠ π/2  (Fig_610.m if branch, p. 432) ──
        sn  = sin(T(π) / n)
        cn  = cos(T(π) / n)
        ctn = cn / sn
        i0  = Complex{T}(kaf)
        is  = im * kaf + (1/n) * ctn + (2/n) * sn / (cn - 1)
        ih  = -im * kaf - (1/n) * ctn + (2/n) * sn / (cn - 1)
        i0, is, ih
    else
        # ── General finite cone  (Fig_610.m else branch, pp. 432–434) ──────────
        # PO contribution = tip integral (from the curved cone surface) + rim edge.
        # PTD adds the fringe-wave edge correction (soft and hard differ at the rim).
        e   = exp(im * 2 * klf)
        tw2 = tan(wf)^2
        sn  = sin(T(π) / n)
        cn  = cos(T(π) / n)
        d1  = cn - one(T)                          # cos(π/n) − 1
        d2  = cn - cos(2 * wf / n)               # cos(π/n) − cos(2w/n)
        tip = im * tw2 * (1 - e) / (2 * kaf)      # PO tip integral
        i0  = tip - tan(wf) * e                   # PO: tip + rim (PO edge)
        is  = tip - (2 * sn / n) * (1/d1 - 1/d2) * e   # PTD soft: tip + PTD edge
        ih  = tip + (2 * sn / n) * (1/d1 + 1/d2) * e   # PTD hard: tip + PTD edge
        i0, is, ih
    end

    return (;
        PO       = abs2(Ish0),   # σ·k²/π, PO approximation
        PTD_soft = abs2(Is),     # σ·k²/π, PTD soft (Dirichlet) BC
        PTD_hard = abs2(Ih),     # σ·k²/π, PTD hard (Neumann) BC
    )
end

# ── Thin wrappers for canonical shapes ────────────────────────────────────────

"""
    paraboloid_rcs(ka, kl, Omega=π/2) -> (; PO, PTD_soft, PTD_hard)

Normalized RCS (σ/(πa²)) for a paraboloid of revolution with rim radius `a` and
axial length `l`, at axial incidence.  Thin wrapper around `scsn_general`.

The rim angle `w = atan(ka / (2·kl))` follows from the paraboloid geometry r² = 2pz
where p = a²/(2l) and the slope at the rim satisfies `tan(w) = p/a = a/(2l)`.
(Ufimtsev App. Ch. 6, `Fig_615.m`, p. 435)

As `kl → 0` (w → π/2) the body degenerates to a flat disk;
as `kl → ∞` (w → 0) it becomes a very elongated nose.  At `kl = 0` exactly,
`atan(Inf) = π/2` is returned, triggering the flat-disk branch automatically.
"""
paraboloid_rcs(ka::Real, kl::Real, Omega::Real = oftype(float(ka), π) / 2) =
    scsn_general(ka, kl, atan(float(ka) / (2 * float(kl))), Omega, :Paraboloid)

"""
    spherical_segment_rcs(ka, w, Omega=π/2) -> (; PO, PTD_soft, PTD_hard)

Normalized RCS (σ/(πa²)) for a spherical-cap segment of a sphere with rim radius `a`
and rim angle `w`.  Thin wrapper around `scsn_general`.

The axial depth `kl = ka·(1−sin(w))/cos(w)` follows from the sphere geometry where
the sphere radius is `R = a/cos(w)`.  (Ufimtsev App. Ch. 6, `Fig_618.m`, p. 436)

- `w → 0`   : deep spherical cap / front hemisphere (σ_PO = πa² exactly for all ka)
- `w → π/2` : shallow cap → flat disk
"""
function spherical_segment_rcs(ka::Real, w::Real, Omega::Real = oftype(float(ka), π) / 2)
    T   = promote_type(typeof(float(ka)), typeof(float(w)), typeof(float(Omega)))
    kaf = T(ka); wf = T(w)
    kl  = kaf * (1 - sin(wf)) / cos(wf)   # axial depth from sphere geometry (Fig_618.m)
    scsn_general(kaf, kl, wf, T(Omega), :Spherical)
end

"""
    flat_disk_rcs(ka, Omega=π/2) -> (; PO, PTD_soft, PTD_hard)

Normalized RCS (σ/(πa²)) for a flat circular PEC disk of radius `a` at normal
incidence.  Thin wrapper: calls `scsn_general` at the flat-disk limit `w = π/2`.

PO result: σ_PO/(πa²) = (ka)² (grows with frequency — large-aperture focusing effect).
(Ufimtsev App. Ch. 6, `Fig_615.m` flat-disk limit, p. 435)
"""
function flat_disk_rcs(ka::Real, Omega::Real = oftype(float(ka), π) / 2)
    T = promote_type(typeof(float(ka)), typeof(float(Omega)))
    scsn_general(T(ka), zero(T), T(π) / 2, T(Omega), :Paraboloid)
end
