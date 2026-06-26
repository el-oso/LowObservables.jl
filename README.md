# LowObservables.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://el-oso.github.io/LowObservables.jl/dev/)
[![CI](https://github.com/el-oso/LowObservables.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/el-oso/LowObservables.jl/actions/workflows/CI.yml)
[![Coverage](https://coveralls.io/repos/github/el-oso/LowObservables.jl/badge.svg?branch=master)](https://coveralls.io/github/el-oso/LowObservables.jl?branch=master)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

**Radar cross-section (RCS) prediction for arbitrary shapes in Julia**, using high-frequency
asymptotics: **Physical Optics (PO)** plus the **Physical Theory of Diffraction (PTD)** edge
correction — the method (originated by Ufimtsev) behind the radar-signature design of the F-117
and B-2. Point it at a surface mesh and get the monostatic RCS.

> ⚠️ Research code under active development. The physics is validated against analytic references
> (see below), but APIs may change.

## What it does

- **Physical Optics solver** — monostatic & bistatic RCS of a perfectly-conducting triangle mesh,
  as a backend-generic [KernelAbstractions](https://github.com/JuliaGPU/KernelAbstractions.jl)
  kernel (CPU now; CUDA backend is a drop-in once GPU hardware is wired up).
- **PTD edge diffraction** — first-order fringe correction from sharp edges (elementary edge waves),
  added coherently to PO. This is what makes RCS accurate at off-specular / edge-on aspects, where PO
  alone collapses into unphysical nulls.
- **Mesh core** — STL/OBJ import, primitive builders (`unit_cube`, `icosphere`, `flat_plate`),
  midpoint refinement, per-face normals/areas/centroids, edge–face adjacency.
- **Exact validation references** — PEC-sphere Mie series and closed-form RCS for canonical bodies
  of revolution (cone, disk, paraboloid, spherical segment), used as the "answer key".
- **Generic over float type** (`Float32`/`Float64`), ready for GPU arrays and automatic differentiation.
- **Interactive documentation** — beginner-friendly, with live WGLMakie widgets (drag the aspect
  angle, watch the RCS and PO-vs-PTD curves update).

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/el-oso/LowObservables.jl")
```

## Quick start

```julia
using LowObservables

λ = 0.03; k = 2π / λ                       # 10 GHz

# --- a flat plate: PO vs PO+PTD ---
plate = flat_plate(0.30, 0.30, 40, 40)    # 30×30 cm, refined
ki = [0.0, 0.0, -1.0]                     # incidence (propagation) direction
ei = [1.0, 0.0, 0.0]                      # E-field polarization ⊥ ki

σ_po  = po_rcs_monostatic(plate, ki, ei; k=k)    # Physical Optics
σ_ptd = ptd_rcs_monostatic(plate, ki, ei; k=k)   # PO + PTD edge correction
to_dbsm(σ_ptd)                                   # → dBsm

# --- load your own shape ---
mesh = load_mesh("aircraft.stl")
ptd_rcs_monostatic(mesh, ki, ei; k=k)

# --- analytic references (validation "answer key") ---
mie_rcs_sphere(0.10, λ) / optical_limit_rcs(0.10)   # PEC sphere → ~1 (optical limit πa²)
```

## Validated physics

Each result below is checked in the test suite against an independent analytic reference:

| Case | Result | Reference |
|------|--------|-----------|
| PEC sphere, high `ka` | σ → π a² | Mie series → optical limit |
| Flat plate, normal incidence | σ → 4π A² / λ² | exact PO of a plate |
| Meshed disk, normal incidence | σ/(πa²) → (ka)² | closed-form `flat_disk_rcs` |
| Flat plate, edge-on | PO ≈ 0, PO+PTD ≠ 0 | rim edge diffraction supplies the RCS |
| Faceted cube, off-specular | sharp 90° edges shift RCS | PTD edge contribution |
| Meshed sphere (smooth) | PO+PTD ≈ PO across refinement | feature-edge gate suppresses facet edges |

## Status

| Phase | | |
|-------|---|---|
| 0 | Mie sphere + interactive-docs infrastructure | ✅ |
| 1 | Closed-form RCS for bodies of revolution | ✅ |
| 2 | Surface-mesh core (import, refine, adjacency) | ✅ |
| 3 | Physical Optics solver (KernelAbstractions) | ✅ |
| 4 | CUDA backend | ⏸ pending GPU hardware |
| 5 | PTD edge diffraction | ✅ |
| 6 | Radar-equation / altitude detectability | ⬜ planned |
| 7 | Shape optimization (RCS minimization) | ⬜ planned |

## Reference

The PO and PTD methods implemented here follow:

> P. Ya. Ufimtsev, *Fundamentals of the Physical Theory of Diffraction*, 2nd ed.,
> Wiley–IEEE Press, 2014. ISBN 978-1-118-75366-8.

Equation and section numbers cited throughout the source code refer to this book.

## License

MIT — see [LICENSE](LICENSE).
