# Scattering Regimes and Bodies of Revolution

This page is for readers who know what RCS is (see [What is RCS?](what-is-rcs.md)) and want
to understand the three scattering regimes, why cones and disks are the canonical test shapes,
and what "axial backscattering" means.

---

## The three scattering regimes

Every RCS problem lives in one of three regimes, determined by the **size parameter** `ka = 2πa/λ`
(wavenumber × characteristic dimension):

| Regime | Condition | Physics | σ behaviour |
|--------|-----------|---------|-------------|
| **Rayleigh** | `ka ≪ 1` | Object is much smaller than a wavelength; acts like an induced dipole | σ ∝ (ka)⁴ — drops steeply with frequency |
| **Resonance** (Mie) | `ka ∼ 1` | Object size comparable to wavelength; standing-wave resonances dominate | σ oscillates above and below πa² |
| **Optical** | `ka ≫ 1` | Short-wave / geometric-optics limit; specular reflection dominates | σ → πa² (geometric cross-section) |

The transition is smooth, but the physics changes completely between regimes.  The Mie series
(Phase 0 of this package) gives an *exact* answer across all three regimes for a sphere.
The formulas on this page operate in the **optical regime** — they use the Physical Theory of
Diffraction (PTD), a high-frequency asymptotic expansion.

---

## Bodies of revolution — why these shapes?

A **body of revolution** (BOR) is a 3-D shape obtained by rotating a 2-D profile around an axis.
Cones, disks, paraboloids, and spherical segments are all BORs.

These shapes are the canonical test cases for three reasons:

1. **Closed-form solutions exist.**  For axial incidence (radar looking straight down the symmetry
   axis), the scattering integral simplifies enough that Ufimtsev's PTD method produces analytical
   formulas.  No numerical integration needed — just evaluate a handful of complex exponentials.

2. **They cover every engineering shape class.**  A cone approximates a missile nose or aircraft
   leading edge.  A flat disk approximates a fuselage panel or inlet aperture.  A paraboloid
   approximates a radome or nose cap.  A spherical segment approximates a curved hull section.

3. **They are the answer key.**  The Phase 1 closed-form solutions here will be used to validate
   the Physical Optics (PO) facet-integration solver coming in Phase 2.  If the numerical PO
   solver agrees with these formulas to within the known PTD correction, the solver is correct.

---

## What does "monostatic axial backscattering" mean?

- **Monostatic**: transmitter and receiver at the same location (standard radar).
- **Axial**: the radar is exactly on the symmetry axis of the body — illuminating the nose head-on.
- **Backscattering**: measuring the power scattered directly back toward the radar.

Axial backscattering is the simplest geometry for a BOR because the problem has rotational symmetry:
the scattered field is the same at every azimuth angle, so only the axial integral survives.
It is also the *worst-case* geometry for most low-observable designs (the nose-on RCS of a
missile or aircraft is usually the most critical).

---

## PO vs PTD: the correction hierarchy

The **Physical Optics (PO)** approximation sums only the contribution from the *illuminated surface*,
treating each patch as a locally flat mirror.  It is simple but ignores **edge diffraction** —
the scattered field from the sharp rim or wedge at the body's edges.

The **Physical Theory of Diffraction (PTD)** adds the fringe-wave correction from the edges:

```
σ_PTD = σ_PO + (edge diffraction term)
```

The correction is different for **soft** (Dirichlet) and **hard** (Neumann) boundary conditions,
which correspond to the two polarisations of the incident wave.

In the closed-form BOR results below, `PTD_soft` and `PTD_hard` straddle the `PO` value — the
soft edge correction typically reduces the backscatter relative to PO, while the hard correction
increases it.  At high ka the corrections oscillate but remain bounded.

---

## Comparison plot: cone axial backscattering  (Ufimtsev Fig. 6.10)

The figure below shows the normalised backscatter cross-section of a 45° half-angle finite cone
(exterior shadow angle Ω = 90°) as a function of the cone-length parameter `kl = ka`.
PO, PTD-soft, and PTD-hard are overlaid.  Directly reproduces the curves in Ufimtsev (2014)
Fig. 6.10, pp. 201–202.

The oscillations come from the interference between the **tip diffraction** and the **base-edge
diffraction**: as `kl` increases, the round-trip phase `2kl` sweeps through multiples of `2π`,
causing constructive and destructive interference.

```@setup cone_plot
using CairoMakie
CairoMakie.activate!()
Makie.inline!(true)
```

```@example cone_plot
using LowObservables, CairoMakie

w     = π / 4      # 45° cone half-angle
Omega = π / 2      # solid cone: shadow boundary at 90°
kls   = range(10.0, 30.0, length = 500)

results  = [cone_rcs(kl, kl, w, Omega) for kl in kls]   # ka = kl for Fig. 6.10
PO_dB    = [10 * log10(r.PO)       for r in results]
soft_dB  = [10 * log10(r.PTD_soft) for r in results]
hard_dB  = [10 * log10(r.PTD_hard) for r in results]

fig = Figure(size = (720, 420))
ax  = Axis(fig[1, 1];
    xlabel = "kl (cone parameter, kl = ka for a right cone)",
    ylabel = "Normalised SCS = σk²/π  (dB)",
    title  = "Finite cone axial backscatter, w = 45°, Ω = 90°  [Ufimtsev Fig. 6.10]",
)
lines!(ax, kls, PO_dB;   color = :black,     linewidth = 1.5, label = "PO")
lines!(ax, kls, soft_dB; color = :steelblue, linewidth = 1.5, linestyle = :dash,
    label = "PTD soft (Dirichlet)")
lines!(ax, kls, hard_dB; color = :tomato,    linewidth = 1.5, linestyle = :dashdot,
    label = "PTD hard (Neumann)")
axislegend(ax; position = :rb)
fig
```

---

## What comes next

Phase 2 will implement a Physical Optics facet-integration solver for arbitrary triangulated
surfaces.  These BOR formulas will serve as the reference:

| Body | `scsn_general` call | Exact PO value (high-ka) |
|------|---------------------|--------------------------|
| Sphere front hemisphere | `scsn_general(ka, ka, 0.0, π/2, :Spherical)` | `PO = 1.0` exactly (σ = πa²) |
| Flat disk | `flat_disk_rcs(ka)` | `PO = (ka)²` (aperture focusing) |
| Paraboloid | `paraboloid_rcs(ka, kl)` | interpolates between above two |

The `PO = 1.0` result for the sphere front hemisphere is an exact analytical identity, not a
numerical approximation — it holds for all `ka`, not just the optical limit.  This is the
primary cross-engine anchor used in the test suite.
