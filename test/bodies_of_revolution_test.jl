# Tests for bodies-of-revolution closed-form RCS formulas.
#
# Reference: Ufimtsev (2014), Appendix to Chapter 6, pp. 431–437.
#
# The most important test is the cross-module anchor: the spherical-segment
# PO result at w=0 must agree with the Mie-series optical limit.

# ── Cross-module anchor: scsn_general (Spherical) ↔ Mie ─────────────────────

@testitem "Cross-module anchor: spherical PO at w=0 equals 1.0 exactly" begin
    # At w = 0 the spherical segment is the front hemisphere of a sphere of radius a
    # (rim at the equator, sphere radius R = a/cos(0) = a, axial depth kl = ka).
    # Analytically:
    #   Ish0 = kR - ka·tan(0)·e^{2ikl} = ka - 0 = ka
    #   PO = |ka|²/ka² = 1.0 EXACTLY for all ka
    #
    # Physical interpretation: σ_PO = πa² for any ka (the PO approximation gives
    # the geometric-optics result exactly for a smooth sphere — no ka dependence).
    using LowObservables
    for ka in [0.5, 1.0, 3π, 10.0, 30.0, 100.0]
        r = scsn_general(ka, ka, 0.0, π / 2, :Spherical)
        @test isapprox(r.PO, 1.0; atol = 1e-12)
        # soft also = 1.0 at w=0 (edge diffraction term cancels: d1 == d2 when 2w=0)
        @test isapprox(r.PTD_soft, 1.0; atol = 1e-12)
    end
end

@testitem "Cross-module anchor: Mie optical limit agrees with scsn_general PO" begin
    # Window-average the Mie series over ka ∈ [30, 31] to suppress Mie ripple.
    # The mean σ/(πa²) should be ≈ 1.0 (optical limit), which matches scsn_general
    # PO = 1.0 at w=0 within 5%.
    #
    # Normalization bridge:
    #   scsn_general returns PO = σ_PO/(πa²)
    #   optical_limit_rcs(a) = πa²
    #   mie_rcs_sphere(a, λ)/optical_limit_rcs(a) = σ_Mie/(πa²)
    # Both should converge to 1.0 as ka → ∞.  This test asserts they agree to 5%.
    using LowObservables
    a   = 1.0
    kas = range(30.0, 31.0, length = 100)
    mie_ratios = [mie_rcs_sphere(a, 2π * a / ka) / optical_limit_rcs(a) for ka in kas]
    mean_mie_ratio = sum(mie_ratios) / length(mie_ratios)

    # scsn_general PO = 1.0 exactly at w=0 for all ka — no ka sweep needed
    scsn_PO = scsn_general(30.0, 30.0, 0.0, π / 2, :Spherical).PO
    @test isapprox(scsn_PO, 1.0; atol = 1e-12)               # exact analytical result
    @test isapprox(mean_mie_ratio, scsn_PO; atol = 0.05)      # Mie → PO within 5%
    @test isapprox(mean_mie_ratio, 1.0; atol = 0.05)          # optical limit check

    # The ratio PO_scsn / (σ_Mie_avg / σ_opt) should be constant = 1.0 ± 5%.
    # This is the cross-engine anchor tying Mie (Phase 0) to BOR formulas (Phase 1).
end

# ── Self-consistency: finite, non-negative, PTD differs from PO ──────────────

@testitem "Self-consistency: PO/PTD finite and non-negative for cone" begin
    # Fig. 6.10 parameter set: w=45°, Omega=90°, ka=kl sweeping 10–30
    using LowObservables
    w     = π / 4
    Omega = π / 2
    for kl in [10.0, 15.0, 20.0, 25.0, 30.0]
        r = cone_rcs(kl, kl, w, Omega)
        @test isfinite(r.PO)       && r.PO       ≥ 0
        @test isfinite(r.PTD_soft) && r.PTD_soft ≥ 0
        @test isfinite(r.PTD_hard) && r.PTD_hard ≥ 0
        # PTD differs from PO (fringe-wave corrections are nonzero for a cone edge)
        @test r.PTD_soft != r.PO
        @test r.PTD_hard != r.PO
    end
end

@testitem "Self-consistency: PO/PTD finite and non-negative for scsn_general" begin
    # Representative parameter sets from Figs. 6.14–6.19
    using LowObservables
    cases = [
        # (ka,   kl,    w,       Omega, shape)
        (3π,  6π,    π/6,   π/2,  :Paraboloid),    # Fig. 6.16 setup
        (3π,  3π,    0.2,   π/2,  :Spherical),      # Fig. 6.18 style
        (3π,  0.0,   π/2,   π/2,  :Paraboloid),     # flat disk limit
        (3π,  3π,    0.0,   π/2,  :Spherical),      # sphere front hemisphere (PO=1)
    ]
    for (ka, kl, w, Omega, shape) in cases
        r = scsn_general(ka, kl, w, Omega, shape)
        @test isfinite(r.PO)       && r.PO       ≥ 0
        @test isfinite(r.PTD_soft) && r.PTD_soft ≥ 0
        @test isfinite(r.PTD_hard) && r.PTD_hard ≥ 0
    end
end

@testitem "PTD fringe corrections are nonzero for hard BC (spherical segment, w=0)" begin
    # At w=0, PO=1 and PTD_soft=1 exactly (tip and rim cancel in the soft case).
    # PTD_hard is non-unity because the hard rim diffraction does NOT cancel.
    using LowObservables
    r = scsn_general(3π, 3π, 0.0, π / 2, :Spherical)
    @test isapprox(r.PO,       1.0; atol = 1e-12)
    @test isapprox(r.PTD_soft, 1.0; atol = 1e-12)
    @test !isapprox(r.PTD_hard, 1.0; atol = 0.01)  # hard BC has nonzero rim diffraction
    @test r.PTD_hard ≥ 0
end

@testitem "Flat disk PO = (ka)²" begin
    # For the flat disk (w=π/2), the PO aperture integral gives σ=4πA²/λ² where A=πa².
    # In the σ/(πa²) normalization: PO = (ka)².
    using LowObservables
    for ka in [1.0, 3.0, 3π, 10.0]
        r = flat_disk_rcs(ka)
        @test isapprox(r.PO, ka^2; rtol = 1e-12)
    end
end

# ── Float32 genericity ────────────────────────────────────────────────────────

@testitem "Float32 in → Float32 out: scsn_general" begin
    using LowObservables
    ka32    = Float32(3π)
    kl32    = Float32(3π)
    w32     = Float32(0.2)
    Omega32 = Float32(π / 2)

    r32 = scsn_general(ka32, kl32, w32, Omega32, :Spherical)
    @test r32.PO       isa Float32
    @test r32.PTD_soft isa Float32
    @test r32.PTD_hard isa Float32

    # Consistent with Float64 to Float32 precision
    r64 = scsn_general(Float64(ka32), Float64(kl32), Float64(w32), Float64(Omega32), :Spherical)
    @test isapprox(Float64(r32.PO),       r64.PO;       rtol = 1e-5)
    @test isapprox(Float64(r32.PTD_soft), r64.PTD_soft; rtol = 1e-5)
    @test isapprox(Float64(r32.PTD_hard), r64.PTD_hard; rtol = 1e-5)
end

@testitem "Float32 in → Float32 out: cone_rcs" begin
    using LowObservables
    r32 = cone_rcs(Float32(10.0), Float32(10.0), Float32(π / 4), Float32(π / 2))
    @test r32.PO       isa Float32
    @test r32.PTD_soft isa Float32
    @test r32.PTD_hard isa Float32

    r64 = cone_rcs(10.0, 10.0, π / 4, π / 2)
    @test isapprox(Float64(r32.PO),       r64.PO;       rtol = 1e-5)
    @test isapprox(Float64(r32.PTD_soft), r64.PTD_soft; rtol = 1e-5)
end

# ── Spot-value sanity checks ─────────────────────────────────────────────────

@testitem "Cone spot values: dB finite, in sane range (Fig. 6.10 params)" begin
    # Fig. 6.10: w=45°, Omega=90°, ka=kl sweeping 10–30.
    # At kl=ka=3π (a multiple of π so exp(2ikl)=1 and the tip integral vanishes):
    #   Ish0 = 0 - tan(π/4)·1 = -1    →  PO = 1   →  0 dB
    # This is an exact analytical value we can pin.
    using LowObservables
    ka  = 3π
    kl  = 3π                        # exp(2i·3π) = exp(6iπ) = 1
    r   = cone_rcs(ka, kl, π / 4, π / 2)
    # PO: tip=0 (since 1 - exp(2ikl) = 0), rim = -tan(π/4) = -1, PO = 1
    @test isapprox(r.PO, 1.0; atol = 1e-10)         # exact: 0 dB
    PO_dB   = 10 * log10(r.PO)
    soft_dB = 10 * log10(r.PTD_soft)
    hard_dB = 10 * log10(r.PTD_hard)
    @test isfinite(PO_dB)   && -20 ≤ PO_dB   ≤ 20
    @test isfinite(soft_dB) && -30 ≤ soft_dB ≤ 20
    @test isfinite(hard_dB) && -10 ≤ hard_dB ≤ 30
    # Sanity: PTD values straddle PO (soft typically lower, hard higher, for a cone)
    # This is a qualitative check — exact values depend on n which depends on w+Omega.
    @test r.PTD_soft < r.PO           # soft diffraction reduces backscatter vs PO
    @test r.PTD_hard > r.PO           # hard diffraction enhances backscatter vs PO
end

@testitem "Cone spot values: representative kl (not multiple of π)" begin
    # At kl=ka=10 (non-trivial phase): just verify the dB value is finite and in [−20, 20].
    # This is a sanity check, not an exact figure pin.
    using LowObservables
    r = cone_rcs(10.0, 10.0, π / 4, π / 2)
    @test isfinite(10 * log10(r.PO))
    @test 0.0 < r.PO < 5.0         # sanity: for a 45° cone, PO ∈ (0, a few) at these params
    @test r.PO ≥ 0 && r.PTD_soft ≥ 0 && r.PTD_hard ≥ 0
end

@testitem "Paraboloid wrapper: consistent with scsn_general" begin
    using LowObservables
    ka = 3π; kl = 6π; Omega = π / 2
    w_expected = atan(ka / (2 * kl))
    r_wrapper = paraboloid_rcs(ka, kl, Omega)
    r_direct  = scsn_general(ka, kl, w_expected, Omega, :Paraboloid)
    @test isapprox(r_wrapper.PO,       r_direct.PO;       rtol = 1e-12)
    @test isapprox(r_wrapper.PTD_soft, r_direct.PTD_soft; rtol = 1e-12)
    @test isapprox(r_wrapper.PTD_hard, r_direct.PTD_hard; rtol = 1e-12)
end

@testitem "Spherical segment wrapper: consistent with scsn_general" begin
    using LowObservables
    ka = 3π; w = 0.3; Omega = π / 2
    kl_expected = ka * (1 - sin(w)) / cos(w)
    r_wrapper = spherical_segment_rcs(ka, w, Omega)
    r_direct  = scsn_general(ka, kl_expected, w, Omega, :Spherical)
    @test isapprox(r_wrapper.PO,       r_direct.PO;       rtol = 1e-12)
    @test isapprox(r_wrapper.PTD_soft, r_direct.PTD_soft; rtol = 1e-12)
    @test isapprox(r_wrapper.PTD_hard, r_direct.PTD_hard; rtol = 1e-12)
end

@testitem "Invalid shape throws ArgumentError" begin
    using LowObservables
    @test_throws ArgumentError scsn_general(3.0, 3.0, 0.2, π / 2, :Cylinder)
end
