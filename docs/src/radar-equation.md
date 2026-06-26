# Radar Range Equation & Stealth

This page explains how far a radar can detect a target, why the detection range
depends on the *fourth root* of RCS, and why making a target stealthy is so hard.

---

## The Radar Range Equation

A monostatic radar transmits power **Pt**, receives the echo, and asks: is there a
target at range **R**?  The answer is in the signal-to-noise ratio (SNR):

```
                Pt · G² · λ² · σ
SNR  =  ─────────────────────────────────────────
         (4π)³ · R⁴ · kB · T₀ · B · F · L
```

| Symbol | Meaning | Typical value |
|--------|---------|---------------|
| Pt     | Transmit power [W] | 1 MW (air-defence radar) |
| G      | Antenna gain (linear) | 10 000 (40 dB) |
| λ = c/f | Wavelength [m] | 3 cm (X-band, 10 GHz) |
| σ      | Target RCS [m²] | 1 m² (fighter-size) |
| R      | Slant range [m] | ? |
| kB     | Boltzmann constant | 1.38 × 10⁻²³ J/K |
| T₀     | Standard noise temperature | 290 K |
| B      | Noise bandwidth [Hz] | 1 MHz |
| F      | Noise figure (linear) | 3 (≈ 4.8 dB) |
| L      | System losses (linear) | 4 (≈ 6 dB) |

The **1/R⁴ factor** is the famous result: two round trips (transmit and receive)
each contribute a 1/R² spreading loss.

---

## Detection Range

Setting SNR = SNR_min and solving for R:

```
R_max = ( (Pt · G² · λ² · σ) / ((4π)³ · kB · T₀ · B · F · L · SNR_min) )^(1/4)
```

For the X-band air-defence radar above with σ = 1 m²:

```julia
julia> using LowObservables
julia> r = RadarSystem()
julia> detection_range(r, 1.0) / 1000   # km
≈ 92 km
```

---

## The σ^(1/4) Problem: Why Stealth is Hard

The detection range scales as **σ^(1/4)** — the *fourth root* of RCS.  This is
the central difficulty for stealth technology.

| RCS reduction | Factor on R_max |
|---------------|-----------------|
| ×10 lower RCS (−10 dB) | ÷ 1.78 |
| ×100 lower RCS (−20 dB) | ÷ 3.16 |
| ×1000 lower RCS (−30 dB) | ÷ 5.62 |
| ×10000 lower RCS (−40 dB) | ÷ 10 |

A 10× RCS reduction cuts detection range by only 1.78×, **not 10×**.
To halve detection range you need a **×16 RCS reduction**.
To cut it by 10× you need a **×10 000 RCS reduction (−40 dB)**.

This is why operational stealth aircraft aim for RCS reductions of 30–40 dB or
more — only then does the range benefit become tactically meaningful.

```julia
julia> detection_range(r, 1.0) / detection_range(r, 0.1)
≈ 1.778   # 10× RCS reduction → only 1.78× range reduction
```

---

## The Radar Horizon

The Earth curves away beneath a target.  Under standard 4/3-earth refraction,
the line-of-sight radar horizon is:

```
R_hor ≈ 4.12 · (√h_target + √h_radar)   [km, h in metres]
```

A target below the horizon is **undetectable regardless of RCS or radar power**.
For a ground-based radar and a target at 10 km altitude:
R_hor ≈ 4.12 × √10 000 = **412 km**.

Low-flying cruise missiles exploit the horizon: flying at 50 m altitude limits
the radar horizon to about 29 km, defeating even a very powerful radar.

---

## Why Stealth Works: Shaping + RAM

### Shaping: redirect the specular flash

For a flat metal plate of area A at normal incidence, Physical Optics gives:

```
σ = (4π / λ²) · A²
```

A 1 m² plate at X-band (λ = 3 cm) has σ = 4π/0.03² ≈ 14 000 m² — the plate
is a **massive** radar reflector when it faces the radar.  Even a small flat facet
on an aircraft contributes a huge RCS spike in the specular direction.

The F-117 Nighthawk solved this by dividing the airframe into flat panels, none of
which faces a typical threat radar.  The specular flash from each panel is
redirected to a small set of known angles — which can be avoided.  The B-2 Spirit
uses a flying-wing shape with no vertical surfaces at all.

### Radar-absorbing materials (RAM)

RAM coatings dissipate radar energy in the material before it can be re-radiated.
RAM is most effective on the leading edges (large RCS contribution from diffraction)
and on surfaces that cannot be shaped away from the radar.

### The Ufimtsev connection

The mathematical foundation for predicting low-observable RCS is the
**Physical Theory of Diffraction (PTD)**, developed by Pyotr Ufimtsev and published
in the Soviet Union in 1962.  The USAF translated it in 1971; Lockheed's Skunk Works
used it to design the F-117, the first operational stealth aircraft.
(Ufimtsev, *Fundamentals of the Physical Theory of Diffraction*, 2nd ed., Introduction
and §1.4.3.)

This package implements PTD edge-wave corrections in `ptd_rcs_monostatic`,
which is used internally by `detection_range_vs_aspect` to compute detection
range as a function of incidence angle.

---

## Interactive: Detection Range vs Aspect Angle

The widget below shows how detection range varies as the aspect angle changes
for a 1 m × 1 m flat PEC plate illuminated by a default X-band air-defence radar.

- **Aspect slider** — rotates the incident direction from broadside (0°, maximum RCS)
  to edge-on (90°, minimum RCS).
- **Stealth-factor slider** — multiplies the base RCS by a constant (emulates RAM or
  shaping that uniformly reduces σ).  Because R_max ∝ σ^(1/4), the curve scales by
  the *fourth root* of the factor.
- **Altitude slider** — sets the target altitude, which determines the radar horizon
  (dashed line).  When the horizon falls below the R_max curve, it becomes the
  binding constraint — extra RCS reduction gives no benefit.

```@setup radar
using WGLMakie, Bonito
WGLMakie.activate!()
Makie.inline!(true)
Page(exportable = true, offline = true)
```

```@example radar
using LowObservables, WGLMakie, Bonito

# ── Precompute ────────────────────────────────────────────────────────────────
# ponytail: precompute base RCS with detection_range_vs_aspect (PTD) once;
#   sliders scale the result analytically via the σ^(1/4) law.

plate = flat_plate(1.0, 1.0, 2, 2)           # 8-face 1 m² PEC plate (z=0 plane)
n_angles = 91
θs   = range(0.0, π / 2; length = n_angles)  # 0° = broadside, 90° = edge-on
dirs = Matrix{Float64}(undef, 3, n_angles)
for (j, θ) in enumerate(θs)
    dirs[:, j] = [-sin(θ), 0.0, -cos(θ)]     # sweep in xz-plane; ki·ẑ<0 for all θ
end
ei    = [0.0, 1.0, 0.0]   # y-polarisation ⊥ ki for all xz-plane incidences ✓
k     = 2π / 0.03          # X-band λ = 3 cm
radar = RadarSystem()      # default X-band air-defence (Pt=1 MW, G=40 dB, f=10 GHz)

base_R_m = detection_range_vs_aspect(plate, dirs, ei; radar, k)   # base R_max [m]

stealth_factors = [1.0, 0.1, 0.01, 0.001]
stealth_labels  = ["1× (no stealth)", "0.1× (−10 dB)", "0.01× (−20 dB)", "0.001× (−30 dB)"]

# GeometryBasics conversion for Makie mesh! recipe
function _to_gb_plate(m)
    pts  = [WGLMakie.Point3f(Float32(m.vertices[1,i]),
                              Float32(m.vertices[2,i]),
                              Float32(m.vertices[3,i])) for i in 1:nvertices(m)]
    tris = [WGLMakie.GLTriangleFace(m.faces[1,j], m.faces[2,j], m.faces[3,j])
            for j in 1:nfaces(m)]
    Makie.GeometryBasics.Mesh(pts, tris)
end
gb_plate = _to_gb_plate(plate)

# ── Widget ────────────────────────────────────────────────────────────────────
App() do session::Session
    sl_aspect  = Bonito.Slider(1:n_angles)   # aspect angle index
    sl_stealth = Bonito.Slider(1:4)          # stealth factor index
    sl_alt     = Bonito.Slider(1:30)         # target altitude [km]

    fig = Figure(size = (680, 760))

    # Top: 3-D view of the flat plate with incident direction arrow
    ax3 = Axis3(fig[1, 1];
        title    = "Plate & incident direction (z = 0 plane)",
        aspect   = :equal,
        azimuth  = 0.4π,
        elevation = 0.25π,
    )
    mesh!(ax3, gb_plate; color = :steelblue, shading = NoShading)

    # Arrow: radar wave arrives from (sin(θ), 0, cos(θ)), points in ki = (-sin, 0, -cos)
    obs_arrow_pt  = map(sl_aspect.value) do j
        [WGLMakie.Point3f(0.5f0 + 1.2f0 * Float32(sin(θs[j])),
                          0.5f0,
                          1.2f0 * Float32(cos(θs[j])))]
    end
    obs_arrow_dir = map(sl_aspect.value) do j
        [WGLMakie.Vec3f(-0.5f0 * Float32(sin(θs[j])),
                         0.0f0,
                        -0.5f0 * Float32(cos(θs[j])))]
    end
    arrows3d!(ax3, obs_arrow_pt, obs_arrow_dir; color = :tomato, shaftradius = 0.04)

    # Bottom: detection range vs aspect angle
    ax2 = Axis(fig[2, 1];
        xlabel = "Aspect angle [°]",
        ylabel = "Detection range [km]",
        title  = "Detection range vs aspect (X-band air-defence)",
    )

    # R_max curve: scales by (stealth_factor)^0.25 when stealth slider moves
    obs_R_pts = map(sl_stealth.value) do si
        sf = stealth_factors[si]^0.25
        [WGLMakie.Point2f(Float32(rad2deg(θs[j])),
                          Float32(base_R_m[j] * sf / 1000))
         for j in 1:n_angles]
    end
    lines!(ax2, obs_R_pts; color = :steelblue, linewidth = 1.5, label = "R_max")

    # Radar horizon: changes with altitude slider
    obs_horizon_pts = map(sl_alt.value) do alt_idx
        h_m = Float64(alt_idx) * 1000.0
        val = Float32(radar_horizon(h_m) / 1000)     # [km]
        [WGLMakie.Point2f(0f0, val), WGLMakie.Point2f(90f0, val)]
    end
    lines!(ax2, obs_horizon_pts;
           color = :tomato, linestyle = :dash, linewidth = 1.5,
           label = "Horizon")
    axislegend(ax2; position = :rt)

    # Moving marker at current aspect + stealth combination
    obs_marker = map(sl_aspect.value, sl_stealth.value) do j, si
        WGLMakie.Point2f(Float32(rad2deg(θs[j])),
                         Float32(base_R_m[j] * stealth_factors[si]^0.25 / 1000))
    end
    scatter!(ax2, obs_marker; color = :black, markersize = 12)

    # Live readout
    obs_label = map(sl_aspect.value, sl_stealth.value, sl_alt.value) do j, si, alt_idx
        R_km  = round(base_R_m[j] * stealth_factors[si]^0.25 / 1000, digits = 1)
        R_hor = round(radar_horizon(alt_idx * 1000.0) / 1000, digits = 1)
        det   = R_km ≤ R_hor ? "DETECTED" : "beyond horizon"
        θ_deg = round(rad2deg(θs[j]), digits = 1)
        "Aspect: $(θ_deg)° | RCS scale: $(stealth_labels[si]) | " *
        "R_max: $(R_km) km | Horizon: $(R_hor) km | $(det)"
    end

    ui = DOM.div(
        "Aspect angle: ", sl_aspect, DOM.br(),
        "Stealth factor (RCS scale): ", sl_stealth, DOM.br(),
        "Target altitude [km]: ", sl_alt, DOM.br(),
        obs_label,
    )
    return Bonito.record_states(session, DOM.div(ui, fig))
end
```
