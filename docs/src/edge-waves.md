# Edge Waves and Fringe Diffraction

This page explains why Physical Optics misses the edge contribution, what the
fringe (nonuniform) current is, and how `fringe_directivity` computes the finite
edge-wave amplitude for any wedge geometry.

**Phase 5b-ii (PO+PTD):** `ptd_rcs_monostatic(mesh, ki, ei; k)` now integrates
the 3-D electromagnetic EEW directivity (`eew_directivity`) along all mesh edges
and adds it coherently to the PO far-field amplitude, so PO+PTD beats PO at
edge-on and near-edge-on aspects where specular return vanishes.

---

## Why PO misses the edge

**Physical Optics (PO)** assumes that the surface current everywhere equals the
*uniform* current induced by the incident plane wave on an infinite flat surface.
This is exact for a smooth, infinite surface but wrong near an edge, where the
current is **concentrated** in a narrow strip and oscillates rapidly.

The missing piece is the **nonuniform (fringe) current** j⁽¹⁾, which is
*localised at the edge* and decays away from it.  PTD adds j⁽¹⁾ to recover the
correct far-field amplitude:

```
E_total = E_PO + E_fringe
```

The fringe current has a *finite* far-field pattern even though the individual
PO and GO terms each diverge at the geometric-optics (GO) boundaries.

---

## The Keller diffraction cone

For a straight edge illuminated at angle γ to the edge direction, the diffracted
rays travel on a **cone** centred on the edge — the **Keller cone** (also called
the diffraction cone or Keller's cone of diffracted rays).  All diffracted rays
make the same angle γ with the edge.  The angular pattern *in the plane perpendicular
to the edge* is given by the 2-D fringe directivity functions f⁽¹⁾ and g⁽¹⁾.

---

## Fringe directivity: f⁽¹⁾ and g⁽¹⁾

For a **perfectly-conducting wedge** of exterior angle α, with incidence angle φ₀
and observation angle φ (both measured in the plane perpendicular to the edge):

- **f⁽¹⁾** (`f1`): soft / E-polarisation (Dirichlet BC, electric field tangent to edge is zero)
- **g⁽¹⁾** (`g1`): hard / H-polarisation (Neumann BC, magnetic field normal to edge is zero)

Both are constructed as `f⁽¹⁾ = f_total − f_PO`, where `f_total` is the exact 2-D
wedge solution and `f_PO` is the PO approximation.  The singularities cancel, leaving
**finite** values at both GO boundaries:

| Boundary | Condition | Physical meaning |
|----------|-----------|-----------------|
| Reflection | φ = π − φ₀ (ψ₂ = π) | Specular reflection from illuminated face |
| Shadow | φ = π + φ₀ (ψ₁ = π) | Edge of the geometric shadow |

Finiteness at these boundaries is the *defining property* of the fringe component —
it cannot be computed as the limit of individual PO terms, which diverge there.

---

## Polar plot: |f⁽¹⁾| and |g⁽¹⁾| for α=270°, φ₀=45°

The figure below shows the fringe directivity magnitudes vs observation angle φ
over the full wedge exterior (0 to 270°).  The two GO boundaries appear as finite
peaks (not poles).  The overall shape reproduces Fig. 4.1 of Ufimtsev (2014).

```@setup fringe_plot
using CairoMakie
CairoMakie.activate!()
Makie.inline!(true)
```

```@example fringe_plot
using LowObservables, CairoMakie

α   = 3π / 2      # 270° exterior wedge angle
φ₀  = π / 4       # 45° incidence (single-side illumination)
φs  = range(1e-3, α - 1e-3, length = 800)

results = [fringe_directivity(φ, φ₀, α) for φ in φs]
f1s     = abs.([r.f1 for r in results])
g1s     = abs.([r.g1 for r in results])

fig = Figure(size = (700, 440))
ax  = Axis(fig[1, 1];
    xlabel = "Observation angle φ [rad]",
    ylabel = "|f⁽¹⁾|, |g⁽¹⁾|",
    title  = "Fringe directivity, α = 270°, φ₀ = 45°  [Ufimtsev Fig. 4.1 shape]",
    xticks = ([0, π/4, π/2, 3π/4, π, 5π/4, 3π/2],
              ["0", "π/4", "π/2", "3π/4", "π", "5π/4", "3π/2"]),
)
lines!(ax, φs, f1s; color = :steelblue, linewidth = 1.5, label = "|f⁽¹⁾| (soft/E-pol)")
lines!(ax, φs, g1s; color = :tomato,    linewidth = 1.5, linestyle = :dash,
    label = "|g⁽¹⁾| (hard/H-pol)")
vlines!(ax, [π - φ₀, π + φ₀]; color = :gray, linewidth = 1, linestyle = :dot)
axislegend(ax; position = :lt)
fig
```

The dotted vertical lines mark the two GO boundaries φ = π − φ₀ and φ = π + φ₀,
where both curves remain **finite** (the PO singularity cancels).

---

---

## 3-D EEW coefficients: F⁽¹⁾ and G⁽¹⁾

The 2-D fringe functions f⁽¹⁾ and g⁽¹⁾ are directivity patterns *in the plane perpendicular
to the edge* — they only tell us the amplitude on the Keller cone.  For computing the RCS of
a real 3-D shape (oblique edges, curved wedges, finite plates), the edge wave must radiate into
**all** observation directions, not just those on the cone.

The electromagnetic **elementary edge wave (EEW)** directivity coefficients F⁽¹⁾ and G⁽¹⁾
generalise f⁽¹⁾ and g⁽¹⁾ to arbitrary polar angle ϑ.  They are vectors in the spherical basis
(ϑ̂, φ̂) with components F_ϑ, F_φ (soft/E-pol) and G_ϑ, G_φ (hard/H-pol).  On the Keller
cone (ϑ = π − γ₀) they reduce to f⁽¹⁾ and g⁽¹⁾ via Eqs. 7.148 and 7.151 of Ufimtsev (2014).

The G_ϑ⁽¹⁾ component is an **exclusively electromagnetic** polarisation-coupling term (Eq. 7.153):
```
G_ϑ(π−γ₀) = [ε(φ₀) − ε(α−φ₀)] · cot(γ₀)
```
where ε(x) = 1 for x ∈ (0,π), else 0.  It is nonzero only for single-face illumination, vanishes
under double-face illumination, and has **no acoustic analog** — it represents coupling between
the two orthogonal polarisations that exists only for the electromagnetic field.

---

## Interactive: PO vs PO+PTD aspect sweep

PO alone misses edge diffraction, so it predicts unphysically low RCS at
off-specular and edge-on aspects; PO+PTD corrects this — the basis of accurate
RCS prediction (and why stealth shaping targets *both* the specular flash and
the edge returns).

The widget below computes monostatic RCS for a 4λ × 4λ flat plate as the aspect
angle sweeps from normal incidence (0°) toward edge-on (90°).  **PO (dashed)**
drops into deep sinc nulls at off-normal aspects; **PO+PTD (solid)** fills them
in with edge-diffraction contributions from the plate rim.

```@setup ptd_compare
using WGLMakie, Bonito
WGLMakie.activate!()
Makie.inline!(true)
Page(exportable = true, offline = true)
```

```@example ptd_compare
using LowObservables, WGLMakie, Bonito

# ── Precompute ────────────────────────────────────────────────────────────────
# ponytail: precompute all 91 angles; only ~160 open rim edges are active in PTD
# (internal flat-plate edges are coplanar → filtered by sharp_threshold).

λ  = 0.1          # wavelength [m]  (X-band, ~10 GHz)
k  = 2π / λ
Lx = 0.4; Ly = 0.4   # 4λ × 4λ plate
mesh = flat_plate(Lx, Ly, 40, 40)   # 3200 faces, λ/10 across
Nf   = nfaces(mesh)

n_angles = 91
θs  = range(0.0, π/2; length = n_angles)
ei  = [0.0, 1.0, 0.0]   # y-pol, always ⊥ ki for xz-plane rotation ✓

σpo_m2  = Vector{Float64}(undef, n_angles)
σptd_m2 = Vector{Float64}(undef, n_angles)
for (j, θ) in enumerate(θs)
    ki = [sin(θ), 0.0, -cos(θ)]
    σpo_m2[j]  = po_rcs_monostatic(mesh, ki, ei; k)
    σptd_m2[j] = ptd_rcs_monostatic(mesh, ki, ei; k)
end
σpo_dBsm  = [to_dbsm(max(s, 1e-15)) for s in σpo_m2]
σptd_dBsm = [to_dbsm(max(s, 1e-15)) for s in σptd_m2]

# Convert plate mesh to GeometryBasics for Makie mesh! recipe
function _plate_to_gb(m)
    pts  = [WGLMakie.Point3f(Float32(m.vertices[1,i]),
                              Float32(m.vertices[2,i]),
                              Float32(m.vertices[3,i])) for i in 1:nvertices(m)]
    tris = [WGLMakie.GLTriangleFace(m.faces[1,j], m.faces[2,j], m.faces[3,j])
            for j in 1:Nf]
    Makie.GeometryBasics.Mesh(pts, tris)
end
gb_plate = _plate_to_gb(mesh)
θ_deg    = rad2deg.(collect(θs))

# ── Widget ────────────────────────────────────────────────────────────────────
App() do session::Session
    sl  = Bonito.Slider(1:n_angles)

    fig = Figure(size = (680, 760))

    # Top: 3-D plate view with incident direction arrow
    ax3 = Axis3(fig[1, 1]; title = "Plate & incident direction",
                aspect = :equal, azimuth = 0.3π, elevation = 0.25π)
    # ponytail: static single-colour mesh. Per-face lit/shadow as Observable
    # Vector{RGBA} is not serializable in offline WGLMakie export.
    mesh!(ax3, gb_plate; color = :steelblue, shading = NoShading)

    # Incident ray arrow: starts in −ki direction, points toward plate centre
    cx = Float32(Lx / 2); cy = Float32(Ly / 2)
    obs_arrow_pt  = map(sl.value) do j
        [WGLMakie.Point3f(cx - 0.8f0*sin(Float32(θs[j])), cy,
                           0.8f0*cos(Float32(θs[j])))]
    end
    obs_arrow_dir = map(sl.value) do j
        [WGLMakie.Vec3f(0.5f0*sin(Float32(θs[j])), 0f0,
                        -0.5f0*cos(Float32(θs[j])))]
    end
    arrows3d!(ax3, obs_arrow_pt, obs_arrow_dir; color = :tomato, shaftradius = 0.015)

    # Bottom: PO vs PO+PTD monostatic RCS
    ax2 = Axis(fig[2, 1];
        xlabel = "Aspect angle θ [°]",
        ylabel = "Monostatic RCS [dBsm]",
        title  = "PO vs PO+PTD — 4λ × 4λ flat plate",
    )
    lines!(ax2, θ_deg, σpo_dBsm;  color = :steelblue, linewidth = 1.5,
           linestyle = :dash,  label = "PO")
    lines!(ax2, θ_deg, σptd_dBsm; color = :tomato,    linewidth = 1.5,
           label = "PO+PTD")
    axislegend(ax2; position = :lb)

    # Moving markers on both curves
    obs_po  = map(sl.value) do j
        WGLMakie.Point2f(Float32(θ_deg[j]), Float32(σpo_dBsm[j]))
    end
    obs_ptd = map(sl.value) do j
        WGLMakie.Point2f(Float32(θ_deg[j]), Float32(σptd_dBsm[j]))
    end
    scatter!(ax2, obs_po;  color = :steelblue, markersize = 12)
    scatter!(ax2, obs_ptd; color = :tomato,    markersize = 12)

    obs_label = map(sl.value) do j
        θd     = round(θ_deg[j], digits = 1)
        po_db  = round(σpo_dBsm[j],  digits = 1)
        ptd_db = round(σptd_dBsm[j], digits = 1)
        "θ = $(θd)°   PO = $(po_db) dBsm   PO+PTD = $(ptd_db) dBsm"
    end
    ui = DOM.div("Aspect angle: ", sl, DOM.br(), obs_label)
    return Bonito.record_states(session, DOM.div(ui, fig))
end
```

**Drag the slider** to sweep the aspect angle from normal incidence (0°) to nearly
edge-on (90°).  Near 0° the two curves coincide — specular scattering dominates
and PTD adds only a small correction.  As θ increases the PO curve falls into
deep sinc nulls; the PO+PTD curve stays elevated because the plate rim acts as a
diffractor and returns energy even when the specular flash has vanished.

---

## API

```@docs
fringe_directivity
eew_directivity
```
