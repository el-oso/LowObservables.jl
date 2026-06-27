# Tutorial: Radar Cross-Section, from Physics to Stealth

This is a guided tour, for readers with first-year calculus and physics, through the ideas
behind `LowObservables.jl`: how radar "sees" an object, why some shapes are nearly invisible
to it, and how to compute and even minimise a radar signature. Each section has a short,
runnable example. By the end you will understand the chain

> **how radar works → how waves scatter → Physical Optics → edge diffraction (PTD) → stealth → shape optimisation.**

If you want to *play* with these ideas, every concept below also has an interactive slider on
its feature page (linked as we go).

```@example tut
using LowObservables
nothing # hide
```

---

## 1. Radar and the radar cross-section

A radar transmits a pulse of radio waves, then listens for the echo. From the echo's delay it
gets range; from its strength it tries to decide *is something there?* The strength of the echo
depends on the target through a single number, the **radar cross-section** (RCS), written ``σ``.

`σ` has units of **area** (m²). The intuition: replace the target with a flat metal disc that
sends the same power back to the radar — the area of that disc is `σ`. It is *not* the physical
size of the object. A stealth jet can have `σ` smaller than a pigeon; a tilted corner of metal
can have `σ` bigger than a truck. Because the numbers span many orders of magnitude, engineers
use a logarithmic unit, **dBsm** (decibels relative to 1 m²): ``σ_\text{dBsm} = 10\log_{10}(σ/1\,\text{m}^2)``.

How does `σ` translate into "can the radar see me?" Through the **radar range equation**. The
echo power falls off as ``1/R^4`` (``1/R^2`` going out, ``1/R^2`` coming back), so the
signal-to-noise ratio is

```math
\text{SNR} \;\propto\; \frac{σ}{R^4}.
```

Turn that around: the maximum **detection range** scales as ``R_\text{max} \propto σ^{1/4}``.
That fourth root is the single most important fact in stealth — to halve the range at which you
are detected, you must cut your RCS by a factor of **16**.

```@example tut
radar = RadarSystem()                       # a generic X-band air-defence radar
R1   = detection_range(radar, 1.0)          # target with σ = 1 m²
R001 = detection_range(radar, 0.01)         # 100× smaller RCS
println("σ = 1.00 m²  → detected at ", round(R1/1e3),   " km")
println("σ = 0.01 m²  → detected at ", round(R001/1e3), " km")
println("100× less RCS bought only ", round(R1/R001, digits=2), "× less range  (≈ 100^¼ = 3.16)")
```

There is also a hard geometric limit: the radar **horizon**. Radio waves travel roughly in
straight lines (bending slightly with the atmosphere), so a low-flying target hides below the
curve of the Earth no matter how large its RCS.

```@example tut
println("A target at 10 km altitude is line-of-sight out to ", round(radar_horizon(10_000.0)/1e3), " km")
```

→ *Explore:* [Radar Equation & Stealth](radar-equation.md) has sliders for altitude and power.

---

## 2. How waves scatter: size versus wavelength

Whether an object scatters strongly, and *how*, is governed by one dimensionless ratio: the
object's size compared with the wavelength ``λ``. It is convenient to use ``ka = 2πa/λ``, where
`a` is a characteristic size — essentially "how many wavelengths across" the object is. Three
regimes appear:

| Regime | Condition | Behaviour |
|--------|-----------|-----------|
| **Rayleigh** | ``ka \ll 1`` (object ≪ wavelength) | ``σ \propto (ka)^4`` — tiny; this is why the sky is blue |
| **Resonance / Mie** | ``ka \sim 1`` | `σ` oscillates as waves wrap around and interfere |
| **Optical** | ``ka \gg 1`` (object ≫ wavelength) | `σ` approaches the geometric value; ray/mirror intuition works |

The perfectly-conducting **sphere** is the one shape we can solve *exactly* (the Mie series,
1908), so it is the benchmark every RCS code is checked against. In the optical regime its RCS
tends to its geometric cross-section, ``σ → πa^2``.

```@example tut
using CairoMakie
a   = 1.0
kas = 10 .^ range(-1, 2, length = 400)         # ka from 0.1 to 100
σn  = [mie_rcs_sphere(a, 2π*a/ka) / (π*a^2) for ka in kas]   # normalised σ/πa²

fig = Figure(size = (672, 360))
ax  = Axis(fig[1,1]; xscale = log10, xlabel = "size parameter  ka",
           ylabel = "σ / πa²", title = "PEC sphere — exact Mie series")
lines!(ax, kas, σn, color = :steelblue)
hlines!(ax, [1.0], color = :tomato, linestyle = :dash, label = "optical limit πa²")
axislegend(ax; position = :rb)
fig
```

Notice the steep Rayleigh rise on the left, the resonance wiggles near ``ka = 1``, and the
flattening toward the optical limit on the right. Radar usually operates in the **optical
regime** (an aircraft is many wavelengths across), which is exactly where the next idea applies.

→ *Explore:* [What is RCS?](what-is-rcs.md) and [Scattering Regimes](scattering-regimes.md).

---

## 3. Physical Optics: scattering as reflection

In the optical regime, waves behave almost like rays, and scattering off a smooth metal surface
looks like **mirror reflection**. Physical Optics (PO) makes this precise. The incident wave
drives electric currents in the metal surface (on the lit side, the induced current is
``\vec J = 2\,\hat n \times \vec H^\text{inc}``); those currents re-radiate, and the scattered
field is the sum of that re-radiation over the lit surface, carefully accounting for the **phase**
each patch contributes.

The phase bookkeeping is the whole game. When a flat panel faces the radar head-on, every patch
re-radiates *in phase* back toward the radar — they add up into an enormous **specular flash**.
Tilt the panel, and the back-toward-radar contributions fall out of phase and largely cancel: the
flash is redirected off to the side. This is the core idea of **shaping**.

```@example tut
λ = 0.03; k = 2π/λ                              # 10 GHz
plate = flat_plate(0.30, 0.30, 30, 30)         # 30 cm × 30 cm metal plate

# head-on
σ_normal = po_rcs_monostatic(plate, [0.0,0,-1.0], [1.0,0,0]; k=k)

# the same plate, but viewed 20° off-normal
θ  = deg2rad(20)
ki = [sin(θ), 0.0, -cos(θ)]; ei = [cos(θ), 0.0, sin(θ)]
σ_tilt = po_rcs_monostatic(plate, ki, ei; k=k)

println("plate head-on : ", round(to_dbsm(σ_normal), digits=1), " dBsm")
println("plate at 20°  : ", round(to_dbsm(σ_tilt),   digits=1), " dBsm")
println("→ a 20° tilt dropped the return by ",
        round(to_dbsm(σ_normal) - to_dbsm(σ_tilt), digits=0), " dB")
```

A 20° tilt cuts the echo by tens of decibels — a factor of thousands. That is why stealth aircraft
are built from flat facets carefully angled *away* from expected radar directions.

(A beautiful side-result of PO: for a large object the *total* scattered power equals twice its
shadow area, ``σ_\text{tot} = 2A`` — the "extinction paradox". You cannot make total scattering
vanish; you can only steer the **backscatter** away from the radar.)

→ *Explore:* [Physical Optics](physical-optics.md) — drag the aspect angle and watch the flash sweep.

---

## 4. Edges diffract: the Physical Theory of Diffraction

PO has a blind spot. It pretends every surface patch is part of an infinite flat mirror, so it
gets the smooth, gently-curved parts right but **mishandles sharp edges**. Yet edges are where
much of a faceted target's return comes from, especially when no facet is angled to give a
specular flash.

What really happens at an edge is **diffraction**: the edge itself acts as a line source, launching
a cone of scattered rays (the *Keller cone*, opening at the same angle the incident ray makes with
the edge). Ufimtsev's **Physical Theory of Diffraction (PTD)** adds exactly the piece PO misses —
the *fringe* field from the edges — so that

```math
\text{(PTD)} \quad σ \;=\; \underbrace{\text{PO}}_{\text{specular}} \;+\; \underbrace{\text{edge diffraction}}_{\text{the fringe correction}}.
```

The difference shows up dramatically off the specular direction. Look at the plate again, almost
**edge-on** (80° off normal), where PO sees almost nothing:

```@example tut
θ  = deg2rad(80)
ki = [sin(θ), 0.0, -cos(θ)]; ei = [cos(θ), 0.0, sin(θ)]
σ_po  = po_rcs_monostatic(plate,  ki, ei; k=k)
σ_ptd = ptd_rcs_monostatic(plate, ki, ei; k=k)
println("edge-on 80°:  PO = ", round(to_dbsm(σ_po),  digits=1),
        " dBsm    PO+PTD = ", round(to_dbsm(σ_ptd), digits=1), " dBsm")
```

PO badly under-predicts the return here; the edge diffraction supplies the missing field. Getting
this right at *every* viewing angle is precisely what made PTD the "industrial-strength" diffraction
theory behind the F-117.

→ *Explore:* [Edge Waves (PTD)](edge-waves.md) — overlay PO and PO+PTD as the aspect sweeps.

---

## 5. Why stealth works

Two levers reduce backscatter:

1. **Shaping** (Section 3): tilt and facet surfaces so the specular flash, and the edge cones,
   point *away* from the threat. This is the dominant effect, and it costs nothing in weight.
2. **Radar-absorbing materials (RAM):** coatings that turn some of the wave's energy into heat,
   reducing the reflected part. RAM helps, but it cannot touch the shadow/forward-scatter
   (the ``σ_\text{tot}=2A`` law), and adds weight and cost.

The F-117 — designed using Ufimtsev's PTD — is the canonical example: a faceted shape that scatters
incoming radar into a few narrow directions that almost never point back at the transmitter, plus
absorbing coatings on the edges.

---

## 6. Designing a low-RCS shape

Instead of shaping by hand, we can let an optimiser do it: deform a mesh to **minimise** its RCS
over a sector of threat directions. There is a catch worth understanding, because it teaches
something general about optimisation.

The cheapest way to lower a plate's RCS is to make the plate *smaller* (``σ \propto A^2``) — a
useless "solution". So the objective must **hold the size fixed** and force the optimiser to
*reshape* rather than shrink. With size pinned, the only way down is to tilt and curve the surface
— exactly the shaping principle, discovered automatically.

```julia
plate = flat_plate(1.0, 1.0, 4, 4)
dirs  = reshape([0.0, 0.0, -1.0], 3, 1)     # threat: head-on
res   = optimize_rcs(plate, dirs, [1.0,0,0]; k = 2π/0.1)
# res.mesh is the reshaped (folded) plate; res.history tracks σ per iteration
```

The optimiser folds the flat plate into a shallow trough that deflects the head-on flash — a huge
PO RCS reduction with the area preserved.

!!! warning "A lesson in trusting your model"
    That folded trough is a *concave corner* — and a corner is a **retroreflector** (it bounces the
    wave back via a double bounce). PO only models a *single* bounce, so it is blind to that return
    and happily reports the trough as nearly invisible. Physically it would be a bright target.
    The optimiser exploited a blind spot in the model. The honest lesson: an optimiser will minimise
    whatever you give it, including the model's errors — real low-observable design needs a scattering
    model with multiple bounces and edge diffraction in the objective, plus constraints against corner
    reflectors. (This is why real stealth shapes are *convex* faceting.)

→ *Explore:* [Shape Optimisation](optimization.md) — scrub the optimiser iterations and watch the fold form.

---

## 7. Under the hood

The methods above are **high-frequency asymptotics**: instead of solving Maxwell's equations
everywhere (impossibly expensive for an object thousands of wavelengths across), they evaluate a
**surface integral** over a triangle mesh of the object — PO over the faces, PTD along the edges.
That surface sum is naturally parallel (each facet is independent), so it runs on a GPU with no
change to the physics — see [GPU (CUDA)](gpu.md). Refining the mesh trades compute for accuracy.

---

## Where to go next

- The interactive feature pages linked throughout let you turn the knobs yourself.
- The definitive reference is **P. Ya. Ufimtsev, *Fundamentals of the Physical Theory of
  Diffraction*, 2nd ed. (Wiley–IEEE, 2014)** — every equation in the source cites it.
- The [API Reference](api.md) lists every function used above.
