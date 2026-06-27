# Stealth shaping: a faceted case study

The defining idea of first-generation stealth (the F-117) is **shaping**: build the aircraft from
flat facets, each angled so that the mirror-like specular flash bounces *away* from where the radar
is — straight ahead. This page demonstrates that with a shape we can actually ship: a synthetic,
all-flat-facet body defined entirely in code, so there are no model licences or restrictions
attached. It is **not** the real F-117 geometry — it is a stylised stand-in that obeys the same
physics.

## The shape

[`faceted_stealth`](@ref) returns a closed, watertight, **F-117-inspired** body: a sharp dart nose,
wings swept ~70° to wingtips at the rear, a flat belly and a ridged faceted top — 9 vertices,
14 facets, every edge shared by exactly two faces. (It is a *stylised stand-in*, not the real
F-117 geometry.)

```@example stealth
using LowObservables

m = faceted_stealth()                       # length=20 m, span=13 m by default
shared = [length(v) for v in values(m.edge_faces)]
println("facets: ", nfaces(m), "   watertight: ", all(==(2), shared),
        "   (Euler V−E+F = ", nvertices(m) - length(m.edge_faces) + nfaces(m), ")")
nothing # hide
```

![faceted stealth dart](assets/stealth_render.png)

Because it is just geometry, you can refine it for finer facets (`refine`) and scale it
(`faceted_stealth(length=…, span=…, height=…)`).

## Its radar signature

We sweep the monostatic RCS around the aircraft in the horizontal plane (azimuth around the vertical
axis, vertical polarisation). At VHF/UHF the facets must be ≲ the wavelength for Physical Optics to be
valid, so we refine to ~3,500 facets and use λ = 1 m — see [Physical Optics](physical-optics.md).

```julia
mr = m; for _ in 1:4; mr = refine(mr); end          # ~3584 facets (≈ 0.5 m, valid at λ=1 m)
λ = 1.0; k = 2π/λ
ϕs = range(0, 2π; length = 181); ei = [0.0, 0.0, 1.0]   # vertical pol ⟂ every horizontal look
σ  = [ptd_rcs_monostatic(mr, [cos(ϕ),sin(ϕ),0], ei; k=k) for ϕ in ϕs]   # ϕ=0 nose-on
```

![faceted stealth RCS](assets/stealth_rcs.png)

| | RCS |
|---|---|
| faceted shape — median | −3.5 dBsm |
| faceted shape — **worst-angle flash** | **+14 dBsm** |
| a flat panel of the same area, at specular | **+53 dBsm** (dashed line) |

## What this shows

- **Faceting kills the specular flash.** A single flat surface facing the radar returns an enormous
  coherent flash — for a panel this size, ``σ = 4πA^2/λ^2 ≈ +53`` dBsm. Break that surface into
  small angled facets and *no* facet presents a large flat area back to the radar, so the worst return
  collapses to **+14 dBsm — about 39 dB (≈8,000×) lower.** That is the core of shaping: never give the
  radar a big mirror.
- **The energy is spread, not removed.** The return becomes a low, lobed pattern (median −3.5 dBsm)
  instead of one giant spike. This is the optical-theorem law ``σ_\text{tot}=2A`` (see
  [Physical Optics](physical-optics.md)): a large object always scatters *something*; shaping spreads
  and redirects that energy rather than eliminating it.
- **Pointing the quiet zones takes optimisation.** This shape is a hand-drawn stand-in, *not* tuned
  for any particular threat direction — its nose-on return is no quieter than its sides. Real stealth
  design carefully angles every facet so the residual lobes miss the expected radar directions; the
  [shape-optimisation](optimization.md) capstone does exactly that automatically, deforming a mesh to
  minimise RCS over a chosen sector.

## Caveats

Treat the **shaping *mechanism*** (no big mirror ⇒ no big flash) as the result, not the precise dBsm:
PO at this facet size is approximate, and a faceted body at VHF is in the resonance regime. The shape
is a stylised, in-code stand-in — for why we don't ship a downloaded F-117 model (licensing/AI-use
restrictions), in-code geometry is simply reproducible and unencumbered.
