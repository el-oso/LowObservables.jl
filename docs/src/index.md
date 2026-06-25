# LowObservables.jl

**Radar cross-section (RCS) simulation for Julia.**

Phase 0 ships the exact Mie series for a perfectly-conducting sphere — the
analytically-exact solution that every RCS simulator is validated against.
Later phases will add physical-optics and geometrical-theory-of-diffraction
facet codes, surface meshes, and GPU acceleration.

## Quick start

```julia
using LowObservables

# Monostatic RCS of a 10 cm metal sphere at 10 GHz (λ ≈ 3 cm)
a = 0.10        # radius [m]
λ = 0.03        # wavelength [m]
σ = mie_rcs_sphere(a, λ)   # → ~0.030 m²

# Sweep over wavelengths
λs = range(0.005, 0.20, length=300)
σs = rcs_vs_size(a, λs)

# Geometric-optics reference
σ_opt = optical_limit_rcs(a)  # = π a²
```

## Model type API

```julia
m = MieSphere(0.10)         # PEC sphere, r = 10 cm
rcs(m, 0.03)                # same as mie_rcs_sphere(0.10, 0.03)
```

## What is RCS?

See [What is RCS?](what-is-rcs.md) for a beginner-friendly explanation
covering units (m²/dBsm), monostatic vs bistatic, and the optical limit.

## API reference

See [API Reference](api.md).
