"""
    LowObservables

Radar cross-section (RCS) simulation for Julia.

Phase 0 delivers the exact Mie series for a perfectly-conducting sphere.
Phase 1 adds closed-form RCS for canonical bodies of revolution (PTD).
Phase 2 adds surface triangle meshes (TriMesh) for arbitrary shapes.
Later phases will add Physical Optics facet codes and GPU acceleration.

## Free sources used
- Bohren & Huffman, *Absorption and Scattering of Light by Small Particles* (1983)
- Wiscombe (1980), Appl. Opt. 19, 1505 — truncation criterion
- Ufimtsev (2014), *Fundamentals of the Physical Theory of Diffraction*, 2nd ed.
- Pharr & Humphreys, *Physically Based Rendering* §3.6 (triangle mesh)
"""
module LowObservables

using SpecialFunctions: besselj, bessely
using TypeContracts
using KernelAbstractions: CPU

include("contracts.jl")
include("mie.jl")
include("bodies_of_revolution.jl")
include("mesh.jl")
include("physical_optics.jl")
include("fringe.jl")
include("elementary_edge_waves.jl")
include("ptd.jl")
include("radar.jl")
include("optimize.jl")

"""
    to_gpu(mesh::TriMesh) -> TriMesh

Move a mesh's per-face arrays (`vertices`, `normals`, `areas`, `centroids`) onto the GPU
so the backend-generic Physical Optics kernel runs on CUDA. Returns a `TriMesh` backed by
`CuArray`s (`faces`/`edge_faces` stay on the CPU — the PO kernel does not read them).

Requires CUDA.jl to be loaded (implemented in `LowObservablesCUDAExt`); calling it without
CUDA loaded raises a `MethodError`.

```julia
using LowObservables, CUDA
σ_gpu = po_rcs_monostatic(to_gpu(mesh), ki, ei; k=k)   # runs on the GPU
```
"""
function to_gpu end

# Phase 0–1: RCS model interface and solvers
export AbstractRCSModel, rcs, rcs_vs_wavelength
export MieSphere
export mie_rcs_sphere, rcs_vs_size, optical_limit_rcs, wiscombe_nmax
export scsn_general, cone_rcs
export paraboloid_rcs, spherical_segment_rcs, flat_disk_rcs

# Phase 2: triangle mesh
export TriMesh
export nvertices, nfaces, total_area
export load_mesh
export unit_cube, icosphere, flat_plate, faceted_stealth
export refine
export plot_mesh

# Phase 3: Physical Optics
export po_rcs, po_rcs_monostatic, po_rcs_sweep
export to_dbsm

# Phase 4: GPU (CUDA backend via extension)
export to_gpu

# Phase 5a: Fringe directivity (PTD edge-diffraction coefficients)
export fringe_directivity

# Phase 5b-i: Electromagnetic EEW directivity (all directions, both polarisations)
export eew_directivity

# Phase 5b-ii: PO + PTD monostatic RCS
export ptd_rcs_monostatic

# Phase 7: Shape optimization
export sector_rcs, optimize_rcs

# Phase 6: Radar range equation / altitude detectability
export RadarSystem
export db2lin, lin2db
export radar_snr, detection_range, radar_horizon, detectable
export detection_range_vs_aspect

end
