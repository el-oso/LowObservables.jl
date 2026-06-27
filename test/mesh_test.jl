# Tests for Phase 2 triangle mesh: TriMesh, primitive builders, refinement,
# adjacency, outward normals, Float32 genericity, and the Makie extension.

# ── Primitive: unit_cube ──────────────────────────────────────────────────────

@testitem "unit_cube: surface area = 6·side²" begin
    using LowObservables

    # side = 1 → total area = 6
    c = unit_cube()
    @test nfaces(c)    == 12
    @test nvertices(c) == 8
    @test isapprox(total_area(c), 6.0; rtol = 1e-14)

    # side = 3 → total area = 54
    c3 = unit_cube(; side = 3.0)
    @test isapprox(total_area(c3), 54.0; rtol = 1e-14)
end

@testitem "unit_cube: face normals are exact axis directions" begin
    using LowObservables, LinearAlgebra
    c = unit_cube()
    # Each face normal should be a unit axis vector (±1 in one coord, 0 elsewhere)
    for i in 1:nfaces(c)
        n = c.normals[:, i]
        @test isapprox(norm(n), 1.0; atol = 1e-14)
        # Exactly one component ≈ ±1, the other two ≈ 0
        extremal = maximum(abs, n)
        @test isapprox(extremal, 1.0; atol = 1e-14)
    end
end

# ── Primitive: icosphere ──────────────────────────────────────────────────────

@testitem "icosphere: face and vertex counts follow 20·4^n, 10·4^n+2" begin
    using LowObservables

    expected_faces    = [20, 80, 320, 1280]
    expected_vertices = [12, 42, 162, 642]
    for (k, n) in enumerate(0:3)
        s = icosphere(n)
        @test nfaces(s)    == expected_faces[k]
        @test nvertices(s) == expected_vertices[k]
    end
end

@testitem "icosphere: area converges toward 4πR² from below" begin
    using LowObservables

    R       = 1.0
    sphere  = 4π * R^2

    areas = [total_area(icosphere(n; R)) for n in 0:3]

    # Areas should be strictly increasing (each subdivision gets closer)
    @test issorted(areas; lt = <)

    # All areas strictly below sphere area
    @test all(a < sphere for a in areas)

    # Level 2 (320 faces) should be within 2 % of 4πR²
    @test isapprox(areas[3], sphere; rtol = 0.02)

    # Level 3 (1280 faces) within 0.5 %
    @test isapprox(areas[4], sphere; rtol = 0.005)
end

@testitem "icosphere: all vertex distances from origin equal R" begin
    using LowObservables, LinearAlgebra
    R = 2.5
    s = icosphere(2; R)
    for i in 1:nvertices(s)
        r = norm(s.vertices[:, i])
        @test isapprox(r, R; rtol = 1e-12)
    end
end

# ── Uniform refinement ────────────────────────────────────────────────────────

@testitem "refine: face count multiplied by exactly 4" begin
    using LowObservables
    for mesh in [unit_cube(), icosphere(1), flat_plate(1.0, 1.0, 3, 3)]
        Nf_before = nfaces(mesh)
        Nf_after  = nfaces(refine(mesh))
        @test Nf_after == 4 * Nf_before
    end
end

@testitem "refine: total area preserved for flat plate" begin
    using LowObservables
    p  = flat_plate(2.0, 3.0, 4, 6)
    r  = refine(p)
    r2 = refine(r)
    @test isapprox(total_area(p),  6.0; rtol = 1e-14)
    @test isapprox(total_area(r),  6.0; rtol = 1e-14)
    @test isapprox(total_area(r2), 6.0; rtol = 1e-14)
end

@testitem "refine: flat subdivision preserves icosphere area" begin
    # `refine` does flat midpoint subdivision (no sphere projection).
    # For an icosphere, each face is an approximately flat planar triangle,
    # so flat subdivision exactly preserves each face's area → total preserved.
    using LowObservables
    s0 = icosphere(0)
    s1 = refine(s0)
    s2 = refine(s1)
    # Area preserved to floating-point precision at each level
    @test isapprox(total_area(s1), total_area(s0); rtol = 1e-12)
    @test isapprox(total_area(s2), total_area(s0); rtol = 1e-10)
    # The icosphere(n) primitive (which DOES project to sphere) converges:
    @test total_area(icosphere(0)) < total_area(icosphere(1)) < total_area(icosphere(2))
    @test isapprox(total_area(icosphere(2)), 4π; rtol = 0.02)
end

@testitem "refine: no duplicate midpoint vertices (shared edges)" begin
    using LowObservables
    # For a cube: 12 edges in base + 12 face diagonals = 18 unique edges.
    # After refine: each edge split → 18 new midpoint verts + 8 original = 26 total.
    # Check with Euler: V - E + F = 2 for a closed orientable surface.
    c  = unit_cube()
    r  = refine(c)
    Nv = nvertices(r)
    Nf = nfaces(r)
    Ne = length(r.edge_faces)    # unique undirected edges
    @test Nv - Ne + Nf == 2      # Euler characteristic
end

# ── Adaptive refinement ───────────────────────────────────────────────────────

@testitem "refine(mask): marked-face count × 4, unmarked kept" begin
    using LowObservables
    p    = flat_plate(1.0, 1.0, 4, 4)
    mask = falses(nfaces(p))
    mask[1:4] .= true         # mark first 4 faces
    r    = refine(p, mask)
    # marked → 4 each (4 × 4 = 16), unmarked kept (28), total = 44
    @test nfaces(r) == nfaces(p) - 4 + 4*4
end

# ── Closed-manifold check ─────────────────────────────────────────────────────

@testitem "unit_cube: every edge shared by exactly 2 faces" begin
    using LowObservables
    c = unit_cube()
    for (_, flist) in c.edge_faces
        @test length(flist) == 2
    end
end

@testitem "icosphere: every edge shared by exactly 2 faces (levels 0-2)" begin
    using LowObservables
    for n in 0:2
        s = icosphere(n)
        for (_, flist) in s.edge_faces
            @test length(flist) == 2
        end
    end
end

# ── Outward normals ───────────────────────────────────────────────────────────

@testitem "unit_cube: signed volume > 0 (outward normals)" begin
    using LowObservables, LinearAlgebra
    c = unit_cube()
    vol = sum(1:nfaces(c)) do i
        i1, i2, i3 = c.faces[:, i]
        v1 = c.vertices[:, i1]; v2 = c.vertices[:, i2]; v3 = c.vertices[:, i3]
        dot(v1, cross(v2, v3))
    end / 6
    @test vol > 0
    @test isapprox(vol, 1.0; rtol = 1e-14)
end

@testitem "icosphere: signed volume > 0 (outward normals)" begin
    using LowObservables, LinearAlgebra
    for n in [0, 1, 2]
        s = icosphere(n)
        vol = sum(1:nfaces(s)) do i
            i1, i2, i3 = s.faces[:, i]
            v1 = s.vertices[:, i1]; v2 = s.vertices[:, i2]; v3 = s.vertices[:, i3]
            dot(v1, cross(v2, v3))
        end / 6
        @test vol > 0
    end
end

@testitem "outward-normal spot check: face normal points away from mesh centroid" begin
    using LowObservables, LinearAlgebra
    for mesh in [unit_cube(), icosphere(1), icosphere(2)]
        mesh_cen = vec(sum(mesh.vertices, dims = 2)) / nvertices(mesh)
        for i in 1:nfaces(mesh)
            face_cen = mesh.centroids[:, i]
            n        = mesh.normals[:, i]
            @test dot(n, face_cen .- mesh_cen) > 0
        end
    end
end

# ── Float32 genericity ────────────────────────────────────────────────────────

@testitem "Float32 mesh: all per-face quantities are Float32" begin
    using LowObservables
    c32  = unit_cube(; T = Float32)
    s32  = icosphere(1; T = Float32)
    p32  = flat_plate(1f0, 1f0, 2, 2; T = Float32)

    for mesh in [c32, s32, p32]
        @test eltype(mesh.vertices)  === Float32
        @test eltype(mesh.normals)   === Float32
        @test eltype(mesh.areas)     === Float32
        @test eltype(mesh.centroids) === Float32
        @test eltype(mesh.faces)     === Int
    end

    # Numerical values should match Float64 to within Float32 precision
    c64  = unit_cube(; T = Float64)
    @test isapprox(Float64(total_area(c32)), total_area(c64); rtol = 1e-6)
end

@testitem "Float32 icosphere area converges same as Float64" begin
    using LowObservables
    sphere32 = Float32(4π)
    a1 = total_area(icosphere(1; T = Float32))
    a2 = total_area(icosphere(2; T = Float32))
    @test a1 < a2 < sphere32
end

# ── Makie extension smoke test ─────────────────────────────────────────────────

@testitem "Makie ext: plot_mesh loads and returns a Figure" tags=[:rendering] begin
    using LowObservables, CairoMakie
    mesh = icosphere(1)
    fig  = plot_mesh(mesh)
    @test fig isa Makie.Figure

    fig2 = plot_mesh(mesh; show_normals = true)
    @test fig2 isa Makie.Figure
end

@testitem "faceted_stealth: watertight F-117-inspired faceted body" begin
    using LowObservables, LinearAlgebra
    m = faceted_stealth()
    @test nvertices(m) == 9
    @test nfaces(m) == 14
    # every undirected edge shared by exactly 2 faces ⇒ closed manifold
    @test all(==(2), [length(v) for v in values(m.edge_faces)])
    @test length(m.edge_faces) == 21          # Euler: V−E+F = 9−21+14 = 2
    @test total_area(m) > 0
    # bounding dimensions honoured exactly (extent-based scaling)
    ext(d) = maximum(m.vertices[d,:]) - minimum(m.vertices[d,:])
    @test isapprox(ext(1), 20.0; rtol=1e-6)   # length
    @test isapprox(ext(2), 13.0; rtol=1e-6)   # span
    @test isapprox(ext(3),  4.0; rtol=1e-6)   # height
    # faces oriented outward: the 5 belly facets all point down
    @test count(<(0), m.normals[3, 1:5]) == 5
    # refinement keeps it watertight (×4 faces, still closed)
    mr = refine(m)
    @test nfaces(mr) == 56
    @test all(==(2), [length(v) for v in values(mr.edge_faces)])
    # Float32 genericity
    m32 = faceted_stealth(T=Float32)
    @test eltype(m32.areas) == Float32
    @test all(==(2), [length(v) for v in values(m32.edge_faces)])
end
