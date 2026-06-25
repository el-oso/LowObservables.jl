# What is Radar Cross-Section?

This page is for complete beginners. If you already know what RCS is, jump
straight to the [API Reference](api.md).

---

## The one-sentence definition

**Radar cross-section (RCS) is the area of a perfect mirror-ball that would
scatter the same amount of power back to the radar as the actual target does.**

It is measured in m², and written σ.

---

## Why "cross-section"?

Think of shining a flashlight at an object. Some of the light bounces
back toward you. If you replaced that object with a perfectly-reflecting
sphere, the sphere that returns the same amount of light has a circular
cross-section of area σ. That is your radar cross-section.

The key insight: σ has nothing to do with the physical size of the object
in any simple way. A stealth aircraft may have an RCS smaller than a bird.
A corner reflector (two planes meeting at 90°) can have an RCS larger than
a truck. Geometry and material matter, not just size.

---

## Units: m² and dBsm

RCS is most naturally measured in **m²** (square metres). In radar
engineering it is often expressed in **dBsm** (decibels relative to one
square metre):

```
σ_dBsm = 10 log₁₀(σ / 1 m²)
```

| σ (m²) | σ (dBsm) | rough equivalent |
|--------|----------|-----------------|
| 0.001  | −30      | insect           |
| 0.01   | −20      | bird             |
| 1      |   0      | small aircraft   |
| 10     | +10      | fighter aircraft |
| 10 000 | +40      | large ship       |

---

## Monostatic vs bistatic

- **Monostatic** (most radars): the transmitter and receiver are at the same
  location (or very close). You measure the power scattered *back* toward the
  radar. This is also called *backscattering*.
- **Bistatic**: transmitter and receiver are at different locations. The
  scattered power in the receiver's direction is measured.

LowObservables.jl currently computes **monostatic** RCS only.

---

## The optical (geometric) limit

For a smooth convex object that is very large compared with the radar
wavelength (so `ka ≫ 1`, where `k = 2π/λ` and `a` is a characteristic
length), the RCS converges to a simple rule:

```
σ → 2 A
```

where `A` is the **silhouette area** — the projected area of the object
as seen from the radar. Factor of 2 comes from the mirror reflection
(a perfectly-reflecting sphere returns twice its physical cross-section).

For a sphere of radius `a`:

```
σ_opt = π a²
```

This is exactly `A = π a²` (the shadow area of a sphere), which equals the
silhouette area. The factor of 2 from reflection and the factor of ½ from
averaging over the hemisphere cancel for a sphere.

You can verify this yourself:

```julia
using LowObservables
a = 1.0   # 1 m sphere
λ = 1e-4  # very short wave (ka ≈ 63 000)
mie_rcs_sphere(a, λ) / optical_limit_rcs(a)  # ≈ 1.0
```

---

## The Rayleigh and resonance regions

| Region          | Condition | Behaviour |
|-----------------|-----------|-----------|
| **Rayleigh**    | `ka ≪ 1`  | σ ∝ (ka)⁴ — extremely small, drops fast |
| **Resonance**   | `ka ~ 1`  | oscillations above *and* below πa² (Mie ripple) |
| **Optical**     | `ka ≫ 1`  | σ → πa² |

The transition is called the **Mie region**, named after Gustav Mie who
worked out the exact series solution in 1908 for a sphere.

---

## The Mie-curve plot (interactive)

The plot below shows the normalised RCS σ/(πa²) as a function of the size
parameter `ka = 2πa/λ`. The horizontal dashed line is the optical limit σ=πa².

**Drag the slider** to move the marker across the curve and watch the
Rayleigh → resonance → optical transition. The widget is a self-contained
WGLMakie figure: the interactivity runs in your browser with no server.

```@setup mie
using WGLMakie, Bonito
WGLMakie.activate!()
Makie.inline!(true)
Page(exportable = true, offline = true)
```

```@example mie
using LowObservables
using WGLMakie, Bonito

a     = 1.0
σ_opt = optical_limit_rcs(a)

# Static curve (dense) + a coarse slider grid the marker snaps to.
kas      = 10 .^ range(-1, 2, length = 400)   # 0.1 … 100 (log)
curve    = [mie_rcs_sphere(a, 2π * a / ka) / σ_opt for ka in kas]
ka_marks = 10 .^ range(-1, 2, length = 60)

App() do session::Session
    sl = Bonito.Slider(1:length(ka_marks))

    fig = Figure(size = (672, 400))
    ax  = Axis(fig[1, 1];
        xscale = log10,
        xlabel = "Size parameter  ka = 2πa/λ",
        ylabel = "Normalised RCS  σ / (πa²)",
        title  = "PEC sphere — Mie series (exact)",
    )
    lines!(ax, kas, curve; color = :steelblue, linewidth = 1.5, label = "Mie")
    hlines!(ax, [1.0]; color = :tomato, linestyle = :dash, linewidth = 1.5,
        label = "Optical limit  σ = πa²")

    # Marker + guide line driven by the slider.
    marker = map(sl.value) do i
        ka = ka_marks[i]
        Point2f(ka, mie_rcs_sphere(a, 2π * a / ka) / σ_opt)
    end
    guide = map(i -> ka_marks[i], sl.value)
    vlines!(ax, guide; color = (:gray, 0.6), linestyle = :dot)
    scatter!(ax, marker; color = :black, markersize = 14)
    axislegend(ax; position = :rb)

    label = map(sl.value) do i
        ka = ka_marks[i]
        r  = mie_rcs_sphere(a, 2π * a / ka) / σ_opt
        "ka = $(round(ka, digits = 2))   σ/(πa²) = $(round(r, digits = 3))"
    end
    ui = DOM.div("Size parameter: ", sl, "  ", label)
    return Bonito.record_states(session, DOM.div(ui, fig))
end
```
