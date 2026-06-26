# Shape Optimisation for Low Observability

This page shows how to deform a mesh to minimise its radar cross-section over a
threat sector — the computational version of F-117-style stealth shaping.

---

## Why shaping works

A flat panel illuminated near-normally creates a **specular flash** toward the
radar.  The RCS of a flat plate at normal incidence is:

```
σ = 4π A² / λ²
```

For a 1 m² plate at X-band (λ = 3 cm), that is **about +52 dBsm** — visible
from hundreds of kilometres.

The fix: **tilt the panel**.  The specular flash obeys the law of reflection.
Tilting by even 10–15° deflects the flash far off the threat axis, reducing
backscatter by many dB.  The F-117 Nighthawk takes this to the extreme:
every surface is a flat panel deliberately angled away from the threat sector,
redirecting all specular returns into a small set of safe directions.

---

## The optimisation problem

Given a mesh with vertex positions **x** (a flat vector of length 3Nv), the
mean monostatic PO RCS over a set of N threat directions is:

```
σ̄(x) = (1/N) Σ_d σ_PO(x, k̂_d)
```

where each `σ_PO` uses the same physical-optics formula as `po_rcs_monostatic`,
but with face geometry (normals, areas, centroids) computed analytically from **x**.

A plain-Julia implementation makes this **ForwardDiff-differentiable**: all face
geometry is smooth in **x**, and the only discontinuity — the lit/shadow indicator
`dot(k̂_i, n̂) < 0` — has zero derivative almost everywhere (the correct
PO shape-derivative).

### The size-preservation constraint

Pure RCS minimisation has a trivial cheat: **shrink the object**.  Since a plate's
σ ∝ A², halving the area cuts σ by 6 dB — without any clever shaping.  A weak
displacement penalty does not stop this; the solver just collapses the mesh.

The fix is to hold the object's **size** with a strong area-preservation term:

```
J(x) = σ̄(x) + μ_area · ((A(x) − A₀)/A₀)²  +  λ_reg · ‖x − x₀‖²
```

where `A(x)` is the total surface area and `A₀` its initial value.  With size
held, shrinking is forbidden and the *only* way to lower σ is to **reshape** —
tilt/curve the surface so the specular flash misses the radar.  The small
`λ_reg · ‖x − x₀‖²` term just anchors the rigid translation modes (which leave σ
unchanged).  For a **closed** body use enclosed volume instead of area, so a
convex shell cannot simply deflate.

---

## Quick start

```julia
using LowObservables

# Build a flat plate (strong specular return at normal incidence)
mesh = flat_plate(1.0, 1.0, 4, 4)        # 32 triangles

# Threat sector: single normal-incidence direction
dirs = reshape([0.0, 0.0, -1.0], 3, 1)   # 3 × 1 matrix
ei   = [1.0, 0.0, 0.0]                   # x-polarisation
k    = 2π / 0.1                           # 10 GHz (λ = 10 cm)

# Initial RCS
x0    = vec(mesh.vertices)
σ_dB  = 10*log10(sector_rcs(x0, mesh.faces, dirs, ei; k))
println("Initial σ = $(round(σ_dB, digits=1)) dBsm")   # ≈ +31 dBsm

# Optimise — size is held by the default area-preservation term (μ_area),
# so the solver must RESHAPE (tilt/curve) rather than shrink.
res  = optimize_rcs(mesh, dirs, ei; k, maxiters = 60)
x_opt = vec(res.mesh.vertices)
σ_opt_dB = 10*log10(sector_rcs(x_opt, res.mesh.faces, dirs, ei; k))
println("Optimised σ = $(round(σ_opt_dB, digits=1)) dBsm")   # ≈ −55 dBsm

# Size preserved, shape changed:
A_ratio = total_area(res.mesh) / total_area(mesh)            # ≈ 1.0 (NOT shrunk)
z_ext   = maximum(res.mesh.vertices[3,:]) - minimum(res.mesh.vertices[3,:])
println("area ratio = $(round(A_ratio, digits=3)),  out-of-plane extent = $(round(z_ext, digits=2)) m")
```

---

## Interactive: watch the mesh reshape

The widget below runs `optimize_rcs` on a 32-face flat plate.  Drag the slider
to scrub through optimiser iterations.  The **wireframe view** shows the flat
plate *folding into a shallow V-trough* — its total area is held fixed, so the
only way to lower σ is to tilt the surface and deflect the specular flash
(red arrow = radar) sideways out of the threat direction.  The **RCS curve**
shows the resulting collapse from +31 dBsm to about −55 dBsm.

```@setup opt
using WGLMakie, Bonito
WGLMakie.activate!()
Makie.inline!(true)
Page(exportable = true, offline = true)
```

```@example opt
using LowObservables, WGLMakie, Bonito

# ── Precompute optimisation history ──────────────────────────────────────────
# ponytail: precompute all states offline; record_states serializes slider states.

mesh0 = flat_plate(1.0, 1.0, 4, 4)   # 32 faces, 25 vertices
dirs  = reshape([0.0, 0.0, -1.0], 3, 1)
ei    = [1.0, 0.0, 0.0]
k     = 2π / 0.1

res    = optimize_rcs(mesh0, dirs, ei; k, maxiters = 60)   # size held by default μ_area
hist_σ = res.history.σ         # Vector{Float64}: sector RCS per iteration [m²]
hist_v = res.history.vertices  # Vector of 3×Nv Float64 vertex matrices
N_iter = length(hist_σ)
F      = mesh0.faces            # face indices fixed across iterations

# Precompute wireframe edge segments per iteration (Observable{Vector{Point3f}})
# ponytail: wireframe via linesegments! — Point3f Observables serialise cleanly
#   in offline Bonito export.  Solid mesh!(obs_mesh) would need a custom WGLMakie
#   binary layout; left for a non-offline (live-server) rendering path.
function _edges(verts, faces)
    pts = WGLMakie.Point3f[]
    for j in axes(faces, 2)
        i1, i2, i3 = faces[1,j], faces[2,j], faces[3,j]
        v1 = WGLMakie.Point3f(Float32(verts[1,i1]), Float32(verts[2,i1]), Float32(verts[3,i1]))
        v2 = WGLMakie.Point3f(Float32(verts[1,i2]), Float32(verts[2,i2]), Float32(verts[3,i2]))
        v3 = WGLMakie.Point3f(Float32(verts[1,i3]), Float32(verts[2,i3]), Float32(verts[3,i3]))
        append!(pts, [v1, v2, v2, v3, v3, v1])
    end
    pts
end
edge_segs = [_edges(hist_v[j], F) for j in 1:N_iter]

σ_dBsm = 10*log10.(max.(hist_σ, 1e-12))   # clamp: σ→0 at the deflecting fold

# ── Widget ────────────────────────────────────────────────────────────────────
App() do session::Session
    sl = Bonito.Slider(1:N_iter)

    fig = Figure(size = (680, 760))

    # Top: 3-D wireframe — updates with slider via Observable{Vector{Point3f}}
    ax3 = Axis3(fig[1, 1]; title = "Mesh deformation (wireframe)",
                aspect = :equal, azimuth = 0.3π, elevation = 0.3π)
    obs_edges = map(sl.value) do j; edge_segs[j]; end
    linesegments!(ax3, obs_edges; color = :steelblue, linewidth = 1.0)

    # Threat direction arrow (radar comes from above: ki = [0,0,-1])
    arrows3d!(ax3,
        [WGLMakie.Point3f(0f0, 0f0, 1.8f0)],
        [WGLMakie.Vec3f( 0f0, 0f0, -0.6f0)];
        color = :tomato, shaftradius = 0.025)

    # Bottom: sector RCS [dBsm] vs iteration
    ax2 = Axis(fig[2, 1];
        xlabel = "Optimiser iteration",
        ylabel = "Sector RCS [dBsm]",
        title  = "RCS reduction (normal incidence, λ = 10 cm)",
    )
    lines!(ax2, 1:N_iter, Float32.(σ_dBsm); color = :steelblue, linewidth = 1.5)
    hlines!(ax2, [Float32(σ_dBsm[end])]; color = :green, linestyle = :dash,
            linewidth = 1, label = "converged")
    axislegend(ax2; position = :rt)

    obs_pt = map(sl.value) do j
        WGLMakie.Point2f(j, Float32(σ_dBsm[j]))
    end
    scatter!(ax2, obs_pt; color = :black, markersize = 12)

    obs_lbl = map(sl.value) do j
        "Iter $(j-1):  σ = $(round(σ_dBsm[j], digits=1)) dBsm"
    end
    ui = DOM.div("Iteration: ", sl, DOM.br(), obs_lbl)
    return Bonito.record_states(session, DOM.div(ui, fig))
end
```

The initial flat plate at normal incidence has σ ≈ +31 dBsm (strong specular return).
After 60 LBFGS iterations the solver reaches σ ≈ −55 dBsm — an **~86 dB reduction**
— by folding the plate into a shallow V-trough that deflects the specular flash off
to the side.  Crucially, the **total area is held within a fraction of a percent of
the original** (A/A₀ ≈ 1.00): this is genuine *reshaping*, not shrinking.

!!! warning "The −55 dBsm is a Physical-Optics number — read this before trusting it"
    The optimiser minimises the **PO** RCS, and PO models only the **single
    (specular) bounce**. The shape it finds is a concave **V-trough**, which deflects
    that single bounce off-axis — so PO reports a near-vanishing return (the 86 dB).

    But a concave dihedral V is a **corner reflector**: it sends the wave straight back
    to the radar via a **double bounce** off its two faces. PO is *blind* to that
    multiple-bounce path, so the V-trough is PO-optimal yet, physically, likely a
    **bright retroreflector** — the opposite of stealth. This is why real low-observable
    shapes (the F-117) use **convex faceting** and scrupulously avoid concave corners
    facing the threat.

    The takeaway is honest and instructive: this demonstrates *gradient-based RCS
    minimisation under a given scattering model*, and that **the optimiser will exploit
    any blindness in that model**. For design use, the objective needs a model that
    includes multiple-bounce returns (and the [PTD edge](edge-waves.md) terms), plus
    constraints that forbid corner-reflector geometry. The machinery here (differentiable
    σ, AD gradients, size-preserved deformation) is exactly what such an extension builds on.

---

## Why not just shrink the shape?

A plate's σ ∝ A², so the cheapest way to lower σ is to make the plate smaller —
the degenerate "collapse to a point" solution.  The size-preservation term
`μ_area · ((A−A₀)/A₀)²` removes that escape: area is pinned, so the optimiser is
forced to *tilt and curve*.  Physical stealth shaping works the same way: an F-117
is roughly the same size as a conventional aircraft — the RCS reduction comes from
orientation, not size.

> **Note — the flat plate is a bifurcation.** A perfectly flat plate at normal
> incidence is symmetric under z → −z (folding up vs down deflects the specular
> equally), so two mirror-image optima exist at identical depth.  Which branch the
> solver picks is set by a tiny symmetry-breaking seed; a slightly non-planar start
> removes the ambiguity.  This is a geometric degeneracy of the demo, not solver
> nondeterminism.

---

## API

| Function | Description |
|----------|-------------|
| `sector_rcs(x, faces, dirs, ei; k)` | Differentiable sector-mean PO RCS from flat vertex vector |
| `optimize_rcs(mesh, dirs, ei; k, λ_reg, maxiters)` | LBFGS shape optimisation; returns `(; mesh, history)` |

See the [API Reference](api.md) for full docstrings.
