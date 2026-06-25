# Tests for Phase 5b-ii: PO + PTD monostatic RCS.
#
# Physics validations:
#   1. Faceted cube: sharp 90° (α=3π/2) edges DIFFRACT — PTD ≠ PO at off-specular aspect
#   2. Icosphere: smooth-body facet edges (α≈π) ARE iterated but suppressed as non-features
#   3. Flat plate, normal incidence: PTD must not corrupt the specular PO peak
#   4. Flat plate, near-edge-on incidence: PTD supplies nonzero diffraction where PO≈0
#   5. Float32 genericity

# ── Test 1: Cube — sharp edges contribute (the key PTD validation) ────────────

@testitem "PTD cube: sharp 90° edges diffract (PTD ≠ PO off-specular)" begin
    using LowObservables
    using LinearAlgebra: normalize

    # Body-diagonal incidence: no face normal is aligned with the radar, so PO has
    # NO specular flash and the 12 sharp 90° edges (α=3π/2) dominate the correction.
    cube = unit_cube(side = 1.0)
    λ = 2.0; k = 2π / λ                       # ka≈1.6: edge diffraction is a sizeable fraction
    ki = -normalize([1.0, 1.0, 1.0])          # propagate toward −[1,1,1]
    ei =  normalize([1.0, -1.0, 0.0])         # ⊥ ki

    σ_po  = po_rcs_monostatic(cube, ki, ei; k)
    σ_ptd = ptd_rcs_monostatic(cube, ki, ei; k)

    @test isfinite(σ_po) && σ_po > 0
    @test isfinite(σ_ptd) && σ_ptd > 0

    # Sharp edges MUST shift the result by a non-trivial margin (was ~12.5%)
    @test abs(σ_ptd - σ_po) / σ_po > 0.05

    # Plausibility: edge correction is a correction, not a blow-up
    # (cube geometric cross-section ~ side²=1; PTD must stay within an order of PO)
    @test 0.3 * σ_po < σ_ptd < 3.0 * σ_po
end

# ── Test 2: Icosphere — facet edges iterated but suppressed as non-features ────

@testitem "PTD icosphere: facet edges iterated, suppressed by feature threshold" begin
    using LowObservables

    a = 1.0
    # icosphere(3): closed, all 2-face edges with α−π ≤ 5.75° (smooth-body tessellation).
    # ka=5: mesh matches analytic PO within 2% (see physical_optics_test.jl).
    sphere = icosphere(3; R = a)
    ki = [0.0, 0.0, -1.0]; ei = [1.0, 0.0, 0.0]
    k  = 5.0

    σ_po  = po_rcs_monostatic(sphere, ki, ei; k)
    σ_ptd = ptd_rcs_monostatic(sphere, ki, ei; k)               # default 20° feature threshold

    @test isfinite(σ_po)  && σ_po  > 0
    @test isfinite(σ_ptd) && σ_ptd > 0

    # With the default feature threshold, the smooth-body facet edges are correctly
    # excluded (not real geometric edges) → PTD ≈ PO.
    @test abs(σ_ptd - σ_po) / σ_po < 0.05

    # Prove the 2-face edges are genuinely ITERATED (not structurally skipped):
    # dropping the threshold to 0 processes every facet edge and DOES change the
    # result substantially.  This is the spurious-faceting sum the threshold removes.
    σ_ptd_all = ptd_rcs_monostatic(sphere, ki, ei; k, sharp_threshold = 0.0)
    @test abs(σ_ptd_all - σ_po) / σ_po > 0.2
end

# ── Test 2: Flat plate, normal incidence — specular peak preserved ─────────────

@testitem "PTD flat plate normal incidence: PTD does not corrupt specular PO" begin
    using LowObservables

    Lx = 1.0; Ly = 1.0
    # 8×8 mesh gives well-resolved facets; use a coarser mesh so test is fast
    plate = flat_plate(Lx, Ly, 8, 8)
    λ = 0.1; k = 2π / λ
    A = total_area(plate)
    σ_exact = 4π * A^2 / λ^2   # exact PO at normal incidence

    ki = [0.0, 0.0, -1.0]; ei = [1.0, 0.0, 0.0]

    σ_po  = po_rcs_monostatic(plate, ki, ei; k)
    σ_ptd = ptd_rcs_monostatic(plate, ki, ei; k)

    # PO alone should be exact at normal incidence (centroid quadrature exact, all phases 0)
    @test isapprox(σ_po, σ_exact; rtol = 1e-10)

    # PO+PTD must not corrupt the peak: within 5% of PO
    @test isfinite(σ_ptd)
    @test isapprox(σ_ptd, σ_po; rtol = 0.05)
end

# ── Test 3: Flat plate, near-edge-on incidence — PTD supplies diffraction ─────

@testitem "PTD flat plate near-edge-on: PTD nonzero where PO≈0" begin
    using LowObservables

    # 85° from normal (5° from edge-on): a fine plate mesh gives PO≈0 (sinc null at
    # grazing), while PTD provides rim half-plane (α=2π) diffraction that stays finite.
    # 16×16 mesh: fine enough that the rapidly-oscillating PO integrand averages to ~0,
    # while the rim-edge PTD sum (only 64 boundary edges) remains physically nonzero.
    plate = flat_plate(1.0, 1.0, 16, 16)
    λ = 0.1; k = 2π / λ

    # ki in xz-plane, 85° from z (normal); ei = y (perpendicular to ki in xz-plane)
    θ = 85.0 * π / 180.0
    ki = [sin(θ), 0.0, -cos(θ)]
    ei = [0.0, 1.0, 0.0]

    σ_po  = po_rcs_monostatic(plate, ki, ei; k)
    σ_ptd = ptd_rcs_monostatic(plate, ki, ei; k)

    # Both must be finite
    @test isfinite(σ_po)
    @test isfinite(σ_ptd)

    # PO must be small (nearly edge-on, sinc-null regime)
    @test σ_po < 0.01      # much less than specular peak ~1257 m²

    # PTD must supply non-trivial rim diffraction and exceed PO
    @test σ_ptd > 1e-2     # finite, nonzero (rim half-plane diffraction ~0.32 m²)
    @test σ_ptd > σ_po     # PTD beats PO at near-edge-on (rim diffraction dominates)
end

# ── Test 4: Float32 genericity ────────────────────────────────────────────────

@testitem "PTD Float32: result is Float32, finite, positive" begin
    using LowObservables

    plate32 = flat_plate(1f0, 1f0, 4, 4; T = Float32)
    ki = Float32[0, 0, -1]
    ei = Float32[1, 0, 0]
    k  = 2f0 * Float32(π) / 0.1f0

    σ = ptd_rcs_monostatic(plate32, ki, ei; k)
    @test σ isa Float32
    @test isfinite(σ)
    @test σ > 0f0

    # Float32 result within 5% of Float64 at normal incidence
    plate64 = flat_plate(1.0, 1.0, 4, 4)
    σ64 = ptd_rcs_monostatic(plate64, [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]; k = Float64(k))
    @test isapprox(Float64(σ), σ64; rtol = 0.05)
end
