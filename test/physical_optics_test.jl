# Tests for Phase 3 Physical Optics solver.
#
# Validation anchors (all @testitem; helpers in @testsetup POHelpers):
#   1. Flat plate, normal incidence: σ = 4πA²/λ²
#   2. Flat plate sinc² null: σ ≈ 0 at first sinc null angle
#   3. Meshed disk vs Phase-1 flat_disk_rcs
#   4. PEC sphere vs analytic PO sphere integral
#   5. Illumination: shadowed incidence → σ = 0
#   6. Float32 end-to-end
#   7. CPU kernel vs reference serial loop (backend genericity)
#   8. bistatic po_rcs and po_rcs_sweep consistency

# ── Shared helpers ────────────────────────────────────────────────────────────

@testsetup module POHelpers

using LowObservables

# Analytic continuous-limit PO RCS for a PEC sphere at monostatic backscatter.
# Derived from the lit-hemisphere integral G_x = 2πa² ∫₀¹ u exp(2ika·u) du.
# Integration by parts: ∫₀¹ u e^{αu} du = ((α−1)e^α + 1)/α²  with α = 2ika.
# σ = (k²/π)|G_x|²  →  σ = (π/(4k²)) |((2ika−1)e^{2ika} + 1)|²
# Large-ka limit: σ → πa²  (PO converges to geometric optics with oscillations).
function po_sphere_analytic(k::Real, a::Real)
    α   = 2 * im * k * a
    val = (α - 1) * exp(α) + 1
    return (π / (4 * k^2)) * abs2(val)
end

# Build a staircase-boundary triangulated disk of radius a by clipping a fine
# rectangular plate to the circle.  Normal incidence PO is exact regardless of
# boundary staircase because all phases vanish at z=0.
function staircase_disk(a::T, n::Int) where {T<:AbstractFloat}
    plate = flat_plate(T(2a), T(2a), n, n)
    mask  = [begin
        cx = plate.centroids[1, i] - a
        cy = plate.centroids[2, i] - a
        cx*cx + cy*cy ≤ a*a
    end for i in 1:nfaces(plate)]
    kept = findall(mask)
    return TriMesh(plate.vertices, plate.faces[:, kept])
end

end  # module POHelpers

# ── Test 1: Plate normal incidence ────────────────────────────────────────────

@testitem "PO plate: normal incidence σ = 4πA²/λ²" begin
    using LowObservables

    λ  = 0.1
    k  = 2π / λ
    plate = flat_plate(1.0, 1.0, 20, 20)
    A = total_area(plate)

    ki = [0.0, 0.0, -1.0]; ei = [1.0, 0.0, 0.0]
    σ       = po_rcs_monostatic(plate, ki, ei; k)
    σ_exact = 4π * A^2 / λ^2

    # centroid quadrature is exact for flat mesh at normal incidence (all phases zero)
    @test isapprox(σ, σ_exact; rtol = 1e-12)

    # to_dbsm conversion
    @test isapprox(to_dbsm(σ), 10 * log10(σ_exact); rtol = 1e-12)
end

# ── Test 2: Plate sinc² null ──────────────────────────────────────────────────

@testitem "PO plate: first sinc null is exactly zero" begin
    using LowObservables

    # Plate in z=0 plane, normal +z. Monostatic sweep in xz-plane.
    # Phase at centroid (cx, cy, 0): k·dot(s−ki, c) = −2k sinθ · cx.
    # Discrete sum vanishes when 2k sinθ · Lx = 2π, i.e. sinθ_null = λ/(2Lx).
    λ  = 0.05; k = 2π / λ
    Lx = 0.5
    plate = flat_plate(Lx, Lx, 30, 30)

    θ_null  = asin(λ / (2Lx))
    ki_null = [sin(θ_null), 0.0, -cos(θ_null)]
    ei_null = [0.0, 1.0, 0.0]   # y-pol ⊥ ki in xz-plane ✓

    σ_null = po_rcs_monostatic(plate, ki_null, ei_null; k)
    @test σ_null < 1e-20    # analytically zero for uniform rectangular mesh

    # Ratio null/peak > 1e15 confirms sinc² pattern
    σ_max = po_rcs_monostatic(plate, [0.0, 0.0, -1.0], ei_null; k)
    @test σ_max / max(σ_null, eps(Float64)) > 1e15
end

# ── Test 3: Disk vs Phase-1 ────────────────────────────────────────────────────

@testitem "PO disk: matches Phase-1 flat_disk_rcs within 5%" setup=[POHelpers] begin
    using LowObservables

    a = 0.5
    λ = 0.1; k = 2π / λ

    disk   = POHelpers.staircase_disk(a, 60)   # 60×60 plate clipped to circle
    A_disk = total_area(disk)

    ki = [0.0, 0.0, -1.0]; ei = [1.0, 0.0, 0.0]
    σ_mesh = po_rcs_monostatic(disk, ki, ei; k)

    # At normal incidence all phases are zero → σ = (k²/π)·A_disk² (exact)
    σ_formula = (k^2 / π) * A_disk^2
    @test isapprox(σ_mesh, σ_formula; rtol = 1e-10)

    # Staircase area approximates πa² to within ~3.5% for a 60×60 grid
    @test isapprox(A_disk, π * a^2; rtol = 0.05)

    # Phase-1 analytic: σ_disk = (ka)² · πa²
    σ_phase1 = flat_disk_rcs(k * a).PO * π * a^2
    @test isapprox(σ_mesh, σ_phase1; rtol = 0.05)
end

# ── Test 4: Sphere vs analytic PO ─────────────────────────────────────────────

@testitem "PO sphere: mesh matches analytic PO integral within 1%" setup=[POHelpers] begin
    using LowObservables

    a = 1.0
    ki = [0.0, 0.0, -1.0]; ei = [1.0, 0.0, 0.0]

    # icosphere(3)=1280 faces at ka=5, icosphere(4)=5120 faces at ka=10
    for (n_sub, ka_val, rtol) in [(3, 5.0, 0.02), (4, 10.0, 0.01)]
        sphere = icosphere(n_sub; R = a)
        k      = ka_val / a
        σ_mesh     = po_rcs_monostatic(sphere, ki, ei; k)
        σ_analytic = POHelpers.po_sphere_analytic(k, a)
        @test isapprox(σ_mesh, σ_analytic; rtol)
    end

    # At ka=10, PO sphere is below πa² (continuous PO oscillates around πa²)
    sphere = icosphere(4; R = a)
    σ_ka10 = po_rcs_monostatic(sphere, ki, ei; k = 10.0)
    @test σ_ka10 < optical_limit_rcs(a)   # 0.91·πa² < πa²

    # Analytic PO converges to πa² at large ka (within 5% at ka=100)
    σ_100 = POHelpers.po_sphere_analytic(100.0, a)
    @test isapprox(σ_100, optical_limit_rcs(a); rtol = 0.05)
end

# ── Test 5a: Illumination — backside shadowed ─────────────────────────────────

@testitem "PO illumination: backside incidence gives σ = 0 exactly" begin
    using LowObservables

    plate = flat_plate(1.0, 1.0, 10, 10)   # normal = +z
    λ = 0.1; k = 2π / λ

    # ki in +z: all faces shadowed (dot([0,0,1],[0,0,1]) = 1 > 0)
    ei = [1.0, 0.0, 0.0]
    @test po_rcs_monostatic(plate, [0.0, 0.0, +1.0], ei; k) == 0.0

    # ki in -z: all faces lit → σ > 0
    @test po_rcs_monostatic(plate, [0.0, 0.0, -1.0], ei; k) > 0.0
end

# ── Test 5b: Illumination — closed convex mesh ────────────────────────────────

@testitem "PO illumination: closed mesh kernel matches manual lit-face filter" begin
    using LowObservables

    sphere = icosphere(2)
    ki = [0.0, 0.0, -1.0]; ei = [1.0, 0.0, 0.0]
    λ = 0.5; k = 2π / λ

    # Kernel result
    σ_full = po_rcs_monostatic(sphere, ki, ei; k)

    # Build mesh with only the lit faces (n_z > 0 for ki=[0,0,-1])
    lit_ids  = findall(i -> sphere.normals[3, i] > 0, 1:nfaces(sphere))
    disk_lit = TriMesh(sphere.vertices, sphere.faces[:, lit_ids])
    σ_lit    = po_rcs_monostatic(disk_lit, ki, ei; k)

    @test 0 < length(lit_ids) < nfaces(sphere)   # some lit, some shadowed
    @test isapprox(σ_full, σ_lit; rtol = 1e-12)  # kernel filters same as explicit mask
end

# ── Test 6: Float32 end-to-end ────────────────────────────────────────────────

@testitem "PO Float32: end-to-end produces Float32 result close to Float64" begin
    using LowObservables

    plate32 = flat_plate(1f0, 1f0, 10, 10; T = Float32)
    ki = Float32[0, 0, -1]; ei = Float32[1, 0, 0]
    k32 = 2f0 * Float32(π) / 0.1f0

    σ32 = po_rcs_monostatic(plate32, ki, ei; k = k32)
    @test σ32 isa Float32
    @test σ32 > 0f0

    # po_rcs_sweep Float32
    dirs = Float32[0 0 0.7071; 0 0 0; -1 -0.7071 -0.7071]
    σs = po_rcs_sweep(plate32, dirs, ei; k = k32)
    @test eltype(σs) === Float32
    @test length(σs) == 3

    # Float32 result within 0.1% of Float64
    plate64 = flat_plate(1.0, 1.0, 10, 10)
    σ64 = po_rcs_monostatic(plate64, [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]; k = Float64(k32))
    @test isapprox(Float64(σ32), σ64; rtol = 1e-3)
end

# ── Test 7: CPU kernel vs reference serial loop ───────────────────────────────

@testitem "PO CPU kernel matches reference serial loop" begin
    using LowObservables

    sphere = icosphere(2)
    λ = 0.3; k = 2π / λ
    ki_raw = [0.5, 0.0, -sqrt(0.75)]; ei_raw = [0.0, 1.0, 0.0]

    σ_kernel = po_rcs_monostatic(sphere, ki_raw, ei_raw; k)

    # Reference: manual serial sum — same arithmetic, no KA kernel
    T = Float64
    _norm3(v) = let n = sqrt(v[1]^2+v[2]^2+v[3]^2); (T(v[1])/n, T(v[2])/n, T(v[3])/n) end
    ki  = _norm3(ki_raw); ei = _norm3(ei_raw)
    s   = (.-ki[1], .-ki[2], .-ki[3])
    kce = (ki[2]*ei[3]-ki[3]*ei[2], ki[3]*ei[1]-ki[1]*ei[3], ki[1]*ei[2]-ki[2]*ei[1])
    smk = (s[1]-ki[1], s[2]-ki[2], s[3]-ki[3])
    Gr1=Gr2=Gr3=Gi1=Gi2=Gi3 = zero(T)
    for i in 1:nfaces(sphere)
        nx=sphere.normals[1,i]; ny=sphere.normals[2,i]; nz=sphere.normals[3,i]
        if ki[1]*nx+ki[2]*ny+ki[3]*nz < zero(T)
            cx=ny*kce[3]-nz*kce[2]; cy=nz*kce[1]-nx*kce[3]; cz=nx*kce[2]-ny*kce[1]
            ph=T(k)*(smk[1]*sphere.centroids[1,i]+smk[2]*sphere.centroids[2,i]+smk[3]*sphere.centroids[3,i])
            sr,si=cos(ph),sin(ph); ai=sphere.areas[i]
            Gr1+=ai*cx*sr; Gi1+=ai*cx*si
            Gr2+=ai*cy*sr; Gi2+=ai*cy*si
            Gr3+=ai*cz*sr; Gi3+=ai*cz*si
        end
    end
    dre=s[1]*Gr1+s[2]*Gr2+s[3]*Gr3; dim=s[1]*Gi1+s[2]*Gi2+s[3]*Gi3
    p1r=Gr1-dre*s[1]; p1i=Gi1-dim*s[1]
    p2r=Gr2-dre*s[2]; p2i=Gi2-dim*s[2]
    p3r=Gr3-dre*s[3]; p3i=Gi3-dim*s[3]
    σ_ref = (T(k)^2/T(π))*(p1r^2+p1i^2+p2r^2+p2i^2+p3r^2+p3i^2)

    @test isapprox(σ_kernel, σ_ref; rtol = 1e-10)
end

# ── Test 8: bistatic and sweep ────────────────────────────────────────────────

@testitem "po_rcs bistatic s=−ki equals po_rcs_monostatic" begin
    using LowObservables

    plate = flat_plate(0.5, 0.5, 8, 8)
    ki = [0.0, 0.0, -1.0]; s = [0.0, 0.0, 1.0]; ei = [1.0, 0.0, 0.0]
    k = 2π / 0.1

    @test isapprox(po_rcs(plate, ki, s, ei; k), po_rcs_monostatic(plate, ki, ei; k); rtol = 1e-14)
end

@testitem "po_rcs_sweep consistent with per-call monostatic" begin
    using LowObservables

    sphere = icosphere(1)
    ei = [0.0, 1.0, 0.0]; k = 2π / 0.2   # y-pol: always ⊥ ki in xz-plane ✓
    θs   = [0.0, π/6, π/3, π/2]
    dirs = hcat([[-sin(θ), 0.0, -cos(θ)] for θ in θs]...)   # 3×4

    σs_sweep = po_rcs_sweep(sphere, dirs, ei; k)
    σs_loop  = [po_rcs_monostatic(sphere, dirs[:, j], ei; k) for j in 1:4]

    @test length(σs_sweep) == 4
    for j in 1:4
        @test isapprox(σs_sweep[j], σs_loop[j]; rtol = 1e-12)
    end
end
