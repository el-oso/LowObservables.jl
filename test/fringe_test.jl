# Tests for fringe directivity functions f⁽¹⁾ and g⁽¹⁾.
#
# Reference: Ufimtsev (2014), §4.3–4.4, pp. 88–100; MATLAB listing fun_fg.m, pp. 447–448.
#
# Key properties verified:
#   1. Finiteness at GO boundaries (the defining property of the fringe wave)
#   2. Eq. 4.29 closed-form values at the reflection boundary
#   3. Float32 genericity

# ── Finiteness at GO boundaries ──────────────────────────────────────────────

@testitem "Fringe: finite at GO boundaries, single-side, α=270°" begin
    using LowObservables
    α  = 3π / 2    # 270° exterior wedge angle, n = 3/2
    φ₀ = π / 4     # 45° incidence — single-side (φ₀ < α−π = π/2)

    # Reflection boundary: φ = π − φ₀, ψ₂ = π
    φ_refl = π - φ₀
    r_refl = fringe_directivity(φ_refl, φ₀, α)
    @test isfinite(r_refl.f1)
    @test isfinite(r_refl.g1)

    # Shadow boundary: φ = π + φ₀, ψ₁ = π
    φ_shad = π + φ₀
    r_shad = fringe_directivity(φ_shad, φ₀, α)
    @test isfinite(r_shad.f1)
    @test isfinite(r_shad.g1)
end

@testitem "Fringe: finite at GO boundaries, single-side, α=300°" begin
    using LowObservables
    α  = 5π / 3    # 300° exterior wedge angle, n = 5/3
    φ₀ = π / 4     # 45° incidence — single-side (φ₀ < α−π = 2π/3)

    φ_refl = π - φ₀
    r_refl = fringe_directivity(φ_refl, φ₀, α)
    @test isfinite(r_refl.f1)
    @test isfinite(r_refl.g1)

    φ_shad = π + φ₀
    r_shad = fringe_directivity(φ_shad, φ₀, α)
    @test isfinite(r_shad.f1)
    @test isfinite(r_shad.g1)
end

@testitem "Fringe: finite at GO boundaries, double-side, α=270°" begin
    # Double-side: φ₀ > α−π = π/2
    using LowObservables
    α  = 3π / 2
    φ₀ = 2π / 3    # 120° incidence — double-side (φ₀ > π/2)

    # ψ₂ = π boundary (reflection, upper face): φ = π − φ₀
    φ_refl = π - φ₀
    r_refl = fringe_directivity(φ_refl, φ₀, α)
    @test isfinite(r_refl.f1)
    @test isfinite(r_refl.g1)

    # 2α−ψ₂ = π boundary (reflection, lower face): ψ₂ = 2α−π → φ = 2α−π−φ₀
    φ_lower = 2α - π - φ₀
    r_lower = fringe_directivity(φ_lower, φ₀, α)
    @test isfinite(r_lower.f1)
    @test isfinite(r_lower.g1)
end

# ── Eq. 4.29 special value at reflection boundary ────────────────────────────

@testitem "Fringe: Eq. 4.29 closed form at φ=π−φ₀ (α=270°, φ₀=45°)" begin
    # At φ = π−φ₀ (ψ₂=π, single-side), Eq. 4.29 gives:
    #   f1 = (1/n)sin(π/n)/(cos(π/n)−cos((π−2φ₀)/n)) + (1/2)cot(φ₀) + (1/2n)cot(π/n)
    #   g1 = (1/n)sin(π/n)/(cos(π/n)−cos((π−2φ₀)/n)) + (1/2)cot(φ₀) − (1/2n)cot(π/n)
    using LowObservables
    α  = 3π / 2;  φ₀ = π / 4;  n = α / π
    φ  = π - φ₀

    common = (1/n) * sin(π/n) / (cos(π/n) - cos((π - 2φ₀) / n)) + 0.5 * cot(φ₀)
    f1_ref = common + (1/(2n)) * cot(π/n)
    g1_ref = common - (1/(2n)) * cot(π/n)

    r = fringe_directivity(φ, φ₀, α)
    @test isapprox(r.f1, f1_ref; rtol = 1e-10)
    @test isapprox(r.g1, g1_ref; rtol = 1e-10)
end

@testitem "Fringe: Eq. 4.29 closed form at φ=π−φ₀ (α=300°, φ₀=30°)" begin
    using LowObservables
    α  = 5π / 3;  φ₀ = π / 6;  n = α / π
    φ  = π - φ₀

    common = (1/n) * sin(π/n) / (cos(π/n) - cos((π - 2φ₀) / n)) + 0.5 * cot(φ₀)
    f1_ref = common + (1/(2n)) * cot(π/n)
    g1_ref = common - (1/(2n)) * cot(π/n)

    r = fringe_directivity(φ, φ₀, α)
    @test isapprox(r.f1, f1_ref; rtol = 1e-10)
    @test isapprox(r.g1, g1_ref; rtol = 1e-10)
end

# ── Float32 genericity ────────────────────────────────────────────────────────

@testitem "Fringe: Float32 in → Float32 out" begin
    using LowObservables
    α32  = Float32(3π / 2)
    φ₀32 = Float32(π / 4)
    φ32  = Float32(π / 3)   # well inside single-side, away from boundaries

    r32 = fringe_directivity(φ32, φ₀32, α32)
    @test r32.f1 isa Float32
    @test r32.g1 isa Float32
    @test isfinite(r32.f1)
    @test isfinite(r32.g1)

    # Consistent with Float64 to within Float32 precision
    r64 = fringe_directivity(Float64(φ32), Float64(φ₀32), Float64(α32))
    @test isapprox(Float64(r32.f1), r64.f1; rtol = 1e-5)
    @test isapprox(Float64(r32.g1), r64.g1; rtol = 1e-5)
end

# ── Self-consistency sweep ────────────────────────────────────────────────────

@testitem "Fringe: finite everywhere in observation sweep, α=270°, φ₀=45°" begin
    # Sweep φ across entire wedge [0, α], including near GO boundaries.
    # All values should be finite (the eps-window branches handle the singularities).
    using LowObservables
    α  = 3π / 2
    φ₀ = π / 4
    φs = range(1e-4, α - 1e-4, length = 1000)
    for φ in φs
        r = fringe_directivity(φ, φ₀, α)
        @test isfinite(r.f1) && isfinite(r.g1)
    end
end
