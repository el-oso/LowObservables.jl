# Physical Optics

This page is for readers new to high-frequency scattering approximations.
If you already know Physical Optics, jump to the [API Reference](api.md).

---

## The key idea: surface currents

When a radar wave hits a metal target, **induced surface currents** flow on the
metal surface and re-radiate energy in all directions.  The exact currents are
hard to compute — they require solving a full electromagnetic boundary-value
problem.

**Physical Optics (PO)** makes one approximation: on the **lit** side of the
target (the half facing the radar), the surface current is the same as it would
be on an infinite flat sheet at that same local orientation.  On the **shadowed**
side, the current is zero.

In equations, the PO surface current at a lit point is:

```
J = 2 n̂ × H_inc
```

where `n̂` is the unit outward normal and `H_inc` is the incident magnetic field.
On the shadow side: `J = 0`.

This is an excellent approximation whenever the wavelength is much shorter than
the local curvature radius of the surface — the **high-frequency** (optical) limit.

---

## Why use a triangle mesh?

The PO integral is a surface integral over the target:

```
G = ∫_lit  (n̂ × (k̂_i × ê_i)) · exp(ik·(ŝ − k̂_i)·r) dS
```

For a real-shaped target (an aircraft, a missile body), this integral has no
closed form.  A **triangle mesh** turns it into a finite sum: one panel per face,
with the centroid approximating the phase integral over that face.

The **centroid approximation** is:

```
∫_face e^{ik(ŝ−k̂_i)·r} dS  ≈  A_i · e^{ik(ŝ−k̂_i)·c_i}
```

where `A_i` is the face area and `c_i` is the centroid.  This is exact when the
phase is constant across the face (e.g., normal incidence on a flat plate), and
converges as the mesh is refined.

---

## Monostatic RCS formula

For **monostatic** (radar transmit = receive), `ŝ = −k̂_i` and the formula
specialises to:

```
σ = (k²/π) |Gperp|²

G      = Σ_lit  A_i (n̂_i × (k̂_i × ê_i)) exp(−2ik k̂_i·c_i)
Gperp  = G − (ŝ·G) ŝ             (component ⊥ to ŝ)
```

with wavenumber `k = 2π/λ`.  This is the standard Physical Optics formula
(Knott, Schaeffer, Tuley, *Radar Cross Section*, §5.2).

---

## Why does a finer mesh give a better answer?

The centroid approximation introduces an error proportional to the **phase
variation** across one face: roughly `k × (face diameter)`.  When this is small
(fine mesh or long wavelength), the approximation is accurate.

`refine(mesh)` splits every triangle into four by inserting edge midpoints.
Four passes multiply the face count by 256 — and the error shrinks by the same
factor.

---

## Where PO fails

Physical Optics assumes the local field equals the flat-sheet solution.  This
breaks down near:

- **Edges** — the current is singular at a sharp edge, not the flat-sheet value.
  The correction is called **Physical Theory of Diffraction (PTD)** — Phase 5.
- **Shadow boundaries** — the transition from lit to dark takes one wavelength,
  not a step function.
- **Concave shapes** — a curved interior can bounce energy back (cavity resonance,
  multiple reflections).  PO only models single-bounce.

---

## Interactive: aspect-angle RCS plot

The widget below shows monostatic RCS (in dBsm) as the radar aspect angle sweeps
from nose-on (0°) to tail-on (180°).  The 3-D view shows which faces are lit
(blue) and which are in shadow (grey) at the current aspect angle.

**Drag the slider** to change the aspect angle.  The shape is an icosphere
(20-face icosahedron, 4 rounds of subdivision = 5120 faces) used as a surrogate
for a rotationally symmetric target.

```@setup po
using WGLMakie, Bonito
WGLMakie.activate!()
Makie.inline!(true)
Page(exportable = true, offline = true)
```

```@example po
using LowObservables, WGLMakie, Bonito

# ── Precompute ────────────────────────────────────────────────────────────────
# ponytail: precompute all 181 angles; record_states captures slider states for
#   offline export.  Live re-run would work but makes the HTML 10× heavier.

a    = 1.0          # sphere radius [m]
λ    = 0.1          # wavelength [m]   (ka = 2πa/λ ≈ 63, well into optics regime)
k    = 2π / λ
mesh = icosphere(4; R = a)   # 5120 faces
Nf   = nfaces(mesh)

# 181 incidence directions sweeping 0°–180° in the xz-plane
n_angles = 181
θs    = range(0.0, π; length = n_angles)
dirs  = Matrix{Float64}(undef, 3, n_angles)
for (j, θ) in enumerate(θs)
    dirs[:, j] = [-sin(θ), 0.0, -cos(θ)]
end
ei = [0.0, 1.0, 0.0]   # y-polarisation ⊥ ki for all xz-plane incidences ✓

σs_m2  = po_rcs_sweep(mesh, dirs, ei; k)          # [m²]
σs_dBsm = to_dbsm.(σs_m2)                         # [dBsm]

# Convert mesh to GeometryBasics for Makie mesh! recipe
function _to_gb(m)
    pts  = [WGLMakie.Point3f(Float32(m.vertices[1,i]),
                              Float32(m.vertices[2,i]),
                              Float32(m.vertices[3,i])) for i in 1:nvertices(m)]
    tris = [WGLMakie.GLTriangleFace(m.faces[1,j], m.faces[2,j], m.faces[3,j])
            for j in 1:Nf]
    Makie.GeometryBasics.Mesh(pts, tris)
end
gb_mesh = _to_gb(mesh)

# ── Widget ────────────────────────────────────────────────────────────────────
App() do session::Session
    sl  = Bonito.Slider(1:n_angles)

    fig = Figure(size = (680, 760))

    # Top: 3-D mesh view (stacked above the RCS plot to fit the docs column width)
    ax3 = Axis3(fig[1, 1]; title = "Shape & incident direction", aspect = :equal,
                azimuth = 0.5π, elevation = 0.2π)
    # ponytail: static single-colour mesh. Per-face lit/shadow recolouring as an
    #   Observable Vector{RGBA} is not serializable as a WGLMakie `vertex_color`
    #   uniform in offline export; a static lit/shadow figure can be added separately.
    mesh!(ax3, gb_mesh; color = :steelblue, shading = NoShading)

    # Incident direction arrow (radar is in the -ki direction)
    obs_arrow_pt  = map(sl.value) do j; [WGLMakie.Point3f( 1.6*sin(θs[j]),  0,  1.6*cos(θs[j]))] end
    obs_arrow_dir = map(sl.value) do j; [WGLMakie.Vec3f(  -0.5*sin(θs[j]),  0, -0.5*cos(θs[j]))] end
    arrows3d!(ax3, obs_arrow_pt, obs_arrow_dir; color = :tomato, shaftradius = 0.02)

    # Bottom: monostatic RCS vs aspect angle
    ax2 = Axis(fig[2, 1];
        xlabel = "Aspect angle [°]",
        ylabel = "Monostatic RCS [dBsm]",
        title  = "PO monostatic sweep",
    )
    lines!(ax2, rad2deg.(θs), σs_dBsm; color = :steelblue, linewidth = 1.5)
    hlines!(ax2, [to_dbsm(optical_limit_rcs(a))]; color = :tomato, linestyle = :dash,
            linewidth = 1, label = "πa² (optical limit)")
    axislegend(ax2; position = :rb)

    # Moving marker on the RCS curve
    obs_marker = map(sl.value) do j
        WGLMakie.Point2f(rad2deg(θs[j]), Float32(σs_dBsm[j]))
    end
    scatter!(ax2, obs_marker; color = :black, markersize = 12)

    # Label
    obs_label = map(sl.value) do j
        θ_deg = round(rad2deg(θs[j]), digits = 1)
        σ_db  = round(σs_dBsm[j], digits = 1)
        "Aspect: $(θ_deg)°   σ = $(σ_db) dBsm"
    end
    ui = DOM.div("Aspect angle: ", sl, DOM.br(), obs_label)
    return Bonito.record_states(session, DOM.div(ui, fig))
end
```

The dashed red line shows the optical-limit reference `σ = πa² ≈ 0 dBsm` for
this 1 m sphere.  The PO result oscillates around this limit as the aspect angle
changes — at broadside (90°) the sphere presents the same circular cross-section
in all directions (for a sphere), while other shapes would show stronger angular
variation.

---

## Quick-start code

```julia
using LowObservables

# Build a mesh
plate = flat_plate(0.5, 0.5, 40, 40)   # 0.5 m × 0.5 m plate, 3200 faces
λ = 0.03                                # 3 cm (≈ 10 GHz, X-band)
k = 2π / λ

# Normal incidence, x-pol
ki = [0.0, 0.0, -1.0]
ei = [1.0, 0.0, 0.0]
σ  = po_rcs_monostatic(plate, ki, ei; k)
println("σ = $(round(σ, digits=2)) m²  ($(round(to_dbsm(σ), digits=1)) dBsm)")

# σ_theory = 4πA²/λ² for a flat plate at normal incidence
A = total_area(plate)
σ_theory = 4π * A^2 / λ^2
println("Theory: $(round(σ_theory, digits=2)) m²")
```

## API

| Function | Description |
|----------|-------------|
| `po_rcs_monostatic(mesh, ki, ei; k)` | Monostatic backscatter RCS [m²] |
| `po_rcs(mesh, ki, s, ei; k)` | Bistatic RCS (arbitrary ŝ) [m²] |
| `po_rcs_sweep(mesh, dirs, ei; k)` | RCS over a matrix of incidence directions |
| `to_dbsm(σ)` | Convert m² to dBsm |

See the [API Reference](api.md) for full docstrings.
