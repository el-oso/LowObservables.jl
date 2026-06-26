# Tests for Phase 6: Radar Range Equation / Altitude Detectability Layer
#
# Physics validations:
#   1. Default X-band radar + σ=1 m² → R_max in [50, 200] km AND matches closed-form
#   2. Scaling laws: ×16 in σ → ×2 in R_max; ×16 in Pt → ×2 in R_max; ×2 in R → /16 SNR
#   3. Stealth 4th-root insight: σ/10 → R_max/10^0.25 ≈ 1.778
#   4. Radar horizon: h_t = 10 km → 412 km
#   5. Detectability gating: beyond R_max OR beyond horizon → false; close high-RCS → true
#   6. Float32 genericity

# ── Test 1: Textbook R_max ────────────────────────────────────────────────────

@testitem "RadarSystem: textbook R_max matches closed-form (σ=1 m²)" begin
    using LowObservables

    r = RadarSystem()   # default X-band air-defence: Pt=1e6 W, G=1e4, f=10 GHz, etc.

    # Independent closed-form recomputation — no magic numbers
    kB = 1.380649e-23; T0 = 290.0; c = 3e8
    λ  = c / r.freq
    R_ref = ((r.Pt * r.G^2 * λ^2 * 1.0) /
             ((4π)^3 * kB * T0 * r.B * r.F * r.L * r.SNR_min))^0.25

    R_lo = detection_range(r, 1.0)

    @test R_lo ≈ R_ref rtol = 1e-10          # matches formula exactly
    @test 50_000 ≤ R_lo ≤ 200_000            # tens to hundreds of km: physically sane
end

# ── Test 2: Scaling laws ──────────────────────────────────────────────────────

@testitem "Scaling laws: σ^(1/4), Pt^(1/4), SNR ∝ 1/R⁴" begin
    using LowObservables

    r = RadarSystem()

    R1  = detection_range(r, 1.0)

    # ×16 in σ → ×2 in R_max  (16^0.25 = 2)
    R16σ = detection_range(r, 16.0)
    @test R16σ ≈ 2 * R1 rtol = 1e-10

    # ×16 in Pt → ×2 in R_max
    r_pt16 = RadarSystem(Pt = 16e6)
    R16Pt  = detection_range(r_pt16, 1.0)
    @test R16Pt ≈ 2 * R1 rtol = 1e-10

    # SNR ∝ 1/R⁴: doubling range → /16 in SNR
    snr1 = radar_snr(r, 1.0, 1_000.0)
    snr2 = radar_snr(r, 1.0, 2_000.0)
    @test snr2 ≈ snr1 / 16 rtol = 1e-10
end

# ── Test 3: Stealth 4th-root insight ──────────────────────────────────────────

@testitem "Stealth insight: σ/10 → R_max / 10^0.25 ≈ 1.778" begin
    using LowObservables

    r   = RadarSystem()
    R1  = detection_range(r, 1.0)
    R01 = detection_range(r, 0.1)   # 10× lower RCS (e.g. RAM treatment)

    # Exact 4th-root: ratio = 10^(1/4) = 10^0.25
    @test R1 / R01 ≈ 10^0.25 rtol = 1e-10
    @test R1 / R01 ≈ 1.77828 rtol = 1e-5   # the famous ~1.78 factor
end

# ── Test 4: Radar horizon ─────────────────────────────────────────────────────

@testitem "Radar horizon: h_target=10 km, h_radar=0 → 412 km" begin
    using LowObservables

    # R_hor = 4.12 * (sqrt(h_t) + sqrt(h_r)) [km], h in [m]
    # = 4120 * sqrt(10000) = 4120 * 100 = 412 000 m
    R_hor = radar_horizon(10_000.0)           # [m]
    @test R_hor ≈ 412_000.0 rtol = 1e-10     # 412 km

    # Two-height formula: radar on a 100 m tower, target at 1 km
    R_hor2 = radar_horizon(1_000.0; h_radar = 100.0)
    expected = 4120.0 * (sqrt(1_000.0) + sqrt(100.0))   # = 4120*(31.62+10) = 171,468 m
    @test R_hor2 ≈ expected rtol = 1e-10
end

# ── Test 5: Detectability gating ──────────────────────────────────────────────

@testitem "Detectability gating: R_max and horizon both gate detection" begin
    using LowObservables

    r = RadarSystem()

    # Case A: high-RCS target, well within R_max and horizon → detectable
    R_maxA = detection_range(r, 100.0)         # large RCS → long R_max
    R_horA = radar_horizon(5_000.0)            # 5 km alt → 4120*70.7 ≈ 291 km
    @test R_maxA > 0 && R_horA > 0
    @test detectable(r, 100.0, 0.5 * R_maxA; h_target = 5_000.0)

    # Case B: beyond R_max (low RCS, long range) → not detectable
    R_maxB = detection_range(r, 0.001)          # stealth target → short R_max
    @test !detectable(r, 0.001, 2.0 * R_maxB; h_target = 5_000.0)

    # Case C: within R_max but beyond horizon (low-altitude target, huge range) → not detectable
    h_low = 100.0                                         # 100 m altitude
    R_hor_low = radar_horizon(h_low)                      # ≈ 4120*10 = 41.2 km
    @test !detectable(r, 100.0, 2.0 * R_hor_low; h_target = h_low)
end

# ── Test 6: Float32 genericity ────────────────────────────────────────────────

@testitem "Float32 genericity: RadarSystem{Float32}, detection_range, radar_snr" begin
    using LowObservables

    r32 = RadarSystem(
        Pt = Float32(1e6), G = Float32(1e4), freq = Float32(10e9),
        B  = Float32(1e6), F = Float32(3.0), L    = Float32(4.0),
        SNR_min = Float32(13.0),
    )
    @test r32 isa RadarSystem{Float32}

    R = detection_range(r32, Float32(1.0))
    @test R isa Float32
    @test R > 0

    snr = radar_snr(r32, Float32(1.0), Float32(1e5))
    @test snr isa Float32
    @test snr > 0

    R_hor = radar_horizon(Float32(10_000.0))
    @test R_hor isa Float32
    @test R_hor ≈ Float32(412_000.0) rtol = 1e-5
end

# ── Test 7: detection_range_vs_aspect smoke test ──────────────────────────────

@testitem "detection_range_vs_aspect: flat plate, 5 aspects, R_max > 0" begin
    using LowObservables
    using LinearAlgebra: normalize

    plate = flat_plate(1.0, 1.0, 2, 2)   # 8-face plate
    λ = 0.03; k = 2π / λ
    θs = range(0.0, π/4; length = 5)
    dirs = hcat([[-sin(θ), 0.0, -cos(θ)] for θ in θs]...)
    ei = [0.0, 1.0, 0.0]

    r = RadarSystem()
    Rs = detection_range_vs_aspect(plate, dirs, ei; radar = r, k)

    @test length(Rs) == 5
    @test all(isfinite, Rs)
    @test all(>(0), Rs)
    # Broadside (θ=0) should give the largest detection range
    @test Rs[1] ≥ Rs[end]
end
