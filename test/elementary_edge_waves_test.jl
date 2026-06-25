# Tests for 3-D electromagnetic EEW directivity coefficients F⁽¹⁾ and G⁽¹⁾.
#
# Reference: Ufimtsev (2014), §7.8.3, Appendix to §7.8.3, MATLAB listing FG.m, pp. 443–445.
#
# Key properties verified:
#   1. On-cone reduction to Phase-5a fringe_directivity (Eqs. 7.148, 7.151) — anchor test
#   2. On-cone polarisation-coupling G_theta (Eq. 7.153) — single-side and double-side
#   3. Float32 genericity (type propagation and Float64 consistency)
#   4. Off-cone smoke test: finite values and F_phi == 0

# ── 1. On-cone anchor: F_theta = −f1/sinγ₀, G_phi = g1/sinγ₀, F_phi = 0 ─────

@testitem "EEW: on-cone reduces to fringe_directivity (Eqs 7.148, 7.151)" begin
    using LowObservables

    configs = [
        # (α,          φ₀,    γ₀,    φ)    — all single-side (φ₀ < α−π)
        (α = 3π/2,    φ₀ = π/4, γ₀ = π/4, φ = π/3),
        (α = 3π/2,    φ₀ = π/4, γ₀ = π/3, φ = π/2),
        (α = 4π/3,    φ₀ = π/6, γ₀ = π/5, φ = π/4),
    ]

    for (; α, φ₀, γ₀, φ) in configs
        ϑ  = π - γ₀                             # exactly on the Keller cone
        eew = eew_directivity(φ, φ₀, γ₀, ϑ, α)
        fg  = fringe_directivity(φ, φ₀, α)

        @test isapprox(eew.F_theta, -fg.f1 / sin(γ₀); rtol = 1e-10)  # Eq.(7.148)
        @test isapprox(eew.G_phi,    fg.g1 / sin(γ₀); rtol = 1e-10)  # Eq.(7.151)
        @test eew.F_phi == 0                                             # Eq.(7.150)
    end
end

# ── 2. On-cone polarisation-coupling G_theta (Eq. 7.153) ─────────────────────

@testitem "EEW: on-cone G_theta = (ε(φ₀)−ε(α−φ₀))·cot(γ₀)  single-side (Eq 7.153)" begin
    using LowObservables

    # Single-side: φ₀ < α−π → ε(φ₀) = 1, ε(α−φ₀) = 0 (α−φ₀ > π) → G_theta = cot(γ₀)
    α = 3π/2;  φ₀ = π/4;  γ₀ = π/4;  φ = π/3
    ϑ = π - γ₀
    eew = eew_directivity(φ, φ₀, γ₀, ϑ, α)

    # ε(π/4)=1, ε(5π/4)=0 → expected = cot(π/4) = 1.0
    @test isapprox(eew.G_theta, cot(γ₀); rtol = 1e-10)

    # Second single-side config with different γ₀
    α2 = 4π/3;  φ₀2 = π/6;  γ₀2 = π/5
    # α2−π = π/3; φ₀2 = π/6 < π/3 → single-side; α2−φ₀2 = 7π/6 > π → ε(α2−φ₀2) = 0
    eew2 = eew_directivity(π/4, φ₀2, γ₀2, π - γ₀2, α2)
    @test isapprox(eew2.G_theta, cot(γ₀2); rtol = 1e-10)
end

@testitem "EEW: on-cone G_theta = 0 for double-side illumination (Eq 7.153)" begin
    using LowObservables

    # Double-side: φ₀ > α−π → ε(φ₀)=1, ε(α−φ₀)=1 → G_theta = 0
    α = 3π/2;  φ₀ = 2π/3;  γ₀ = π/4;  φ = π/2  # φ₀=120° > α−π=90° → double-side
    # α−φ₀ = 5π/6 ∈ (0,π) → ε(α−φ₀)=1; ε(φ₀)=1 → diff = 0
    eew = eew_directivity(φ, φ₀, γ₀, π - γ₀, α)
    @test eew.G_theta == 0
end

# ── 3. Float32 genericity ────────────────────────────────────────────────────

@testitem "EEW: Float32 in → Float32 out, on-cone and off-cone" begin
    using LowObservables

    α32  = Float32(3π/2)
    φ₀32 = Float32(π/4)
    γ₀32 = Float32(π/4)
    φ32  = Float32(π/3)

    # On-cone
    ϑ_on = Float32(π) - γ₀32   # triggers on-cone branch exactly
    r_on = eew_directivity(φ32, φ₀32, γ₀32, ϑ_on, α32)
    @test r_on.F_theta isa Float32
    @test r_on.F_phi   isa Float32
    @test r_on.G_theta isa Float32
    @test r_on.G_phi   isa Float32
    @test isfinite(r_on.F_theta) && isfinite(r_on.G_phi)

    # Consistent with Float64 on-cone formula to within Float32 precision.
    # NOTE: Float64(Float32(π)-Float32(γ₀)) ≠ Float64(π)-Float64(γ₀) by ~3e-8,
    # so the Float64 eew_directivity would take the off-cone branch. Use fringe_directivity
    # directly (the on-cone Eqs. 7.148/7.151) as the reference instead.
    fg64  = fringe_directivity(Float64(φ32), Float64(φ₀32), Float64(α32))
    γ₀64  = Float64(γ₀32)
    @test isapprox(Float64(r_on.F_theta), -fg64.f1 / sin(γ₀64); rtol = 1e-5)
    @test isapprox(Float64(r_on.G_phi),    fg64.g1 / sin(γ₀64); rtol = 1e-5)

    # Off-cone (ϑ = π/2, cone is at 3π/4). Off the cone the EEW coefficients are
    # COMPLEX in general (sigma12 Eq.7.77), so components are ComplexF32, not Float32.
    ϑ_off = Float32(π/2)
    r_off = eew_directivity(φ32, φ₀32, γ₀32, ϑ_off, α32)
    @test r_off.F_theta isa ComplexF32
    @test r_off.G_theta isa ComplexF32
    @test r_off.G_phi   isa ComplexF32
    @test isfinite(r_off.F_theta) && isfinite(r_off.G_theta) && isfinite(r_off.G_phi)
    @test r_off.F_phi == 0
end

# ── 5. Continuity across the Keller cone (exercises the off-cone U/V branch) ──

@testitem "EEW: continuous across the Keller cone (off-cone U/V branch)" begin
    using LowObservables

    # NOTE: φ is kept away from α/2 so that β₂ stays clear of the sigma12 fold βₖ = 2γ₀
    # (Eqs 7.76↔7.77). On that fold the published sigma12 switches its real↔complex branch
    # and is itself discontinuous — a known artifact of the book's algorithm, unrelated to
    # the cone reduction these configs verify.
    configs = [
        (α = 3π/2,  φ₀ = π/4, γ₀ = π/4, φ = π/3),
        (α = 3π/2,  φ₀ = π/3, γ₀ = π/3, φ = 2π/3),
        (α = 4π/3,  φ₀ = π/5, γ₀ = π/4, φ = 2π/3),
    ]

    for (; α, φ₀, γ₀, φ) in configs
        on = eew_directivity(φ, φ₀, γ₀, π - γ₀, α)            # on-cone (calls fun_fg)
        for dϑ in (1e-4, -1e-4)
            off = eew_directivity(φ, φ₀, γ₀, π - γ₀ + dϑ, α)  # off-cone (U/V branch)
            # Continuous function ⇒ |off − on| ~ O(dϑ); a broken branch jumps by O(1).
            @test isapprox(off.F_theta, on.F_theta; atol = 1e-2)
            @test isapprox(off.G_phi,   on.G_phi;   atol = 1e-2)
            @test isapprox(off.G_theta, on.G_theta; atol = 1e-2)
            @test off.F_phi == 0
        end
    end
end

# ── 4. Off-cone smoke test: finite values, F_phi = 0 ─────────────────────────

@testitem "EEW: off-cone returns finite values and F_phi = 0" begin
    using LowObservables

    # Well away from singularities and GO boundaries
    α  = 3π/2;  φ₀ = π/4;  γ₀ = π/4;  φ = π/3

    for ϑ in [π/6, π/3, π/2, 2π/3, 5π/6]   # various off-cone observation angles
        eew = eew_directivity(φ, φ₀, γ₀, ϑ, α)
        @test isfinite(eew.F_theta)
        @test isfinite(eew.G_theta)
        @test isfinite(eew.G_phi)
        @test eew.F_phi == 0
    end
end
