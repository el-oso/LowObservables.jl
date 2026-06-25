"""
Surface triangle mesh for Physical Optics (Phase 2).

`TriMesh{T,A}` is a struct-of-arrays mesh parametric over the float type `T`
and the float matrix type `A<:AbstractMatrix{T}`.  On CPU, `A = Matrix{T}`.
Phase 3 GPU kernels (KernelAbstractions) promote `A` to `CuMatrix{T}` by
wrapping the CPU construction output — see the ponytail note on the type.

References:
- Pharr & Humphreys, *Physically Based Rendering* §3.6 (triangle mesh)
- Koenderink (1990), solid-angle / divergence-theorem signed-volume formula
"""

using LinearAlgebra: cross, dot, norm
using StaticArrays: SVector
using GeometryBasics
using FileIO

# ── TriMesh type ──────────────────────────────────────────────────────────────

"""
    TriMesh{T, A<:AbstractMatrix{T}}

Surface triangle mesh as struct-of-arrays.

# Fields (column-major)
- `vertices   :: A`           — 3×Nv; column k = k-th vertex position [m]
- `faces      :: Matrix{Int}` — 3×Nf; column k = 1-based vertex indices
- `normals    :: A`           — 3×Nf; unit outward face normals
- `areas      :: Vector{T}`   — Nf; face areas [m²]
- `centroids  :: A`           — 3×Nf; face centroids [m]
- `edge_faces :: Dict{Tuple{Int,Int},Vector{Int}}` — undirected edge → face indices

# ponytail: `faces` (Matrix{Int}) and `areas` (Vector{T}) are CPU-only types.
# Phase 3 GPU path adds an integer-array type parameter for `faces` (CuMatrix{Int32})
# and promotes `areas` to a matching 1D CuArray.  Until then, construct on CPU via
# `TriMesh(verts, faces)` and move floats to GPU: `CuMatrix(mesh.vertices)`, etc.

# ponytail: `edge_faces` is CPU Dict preprocessing for Phase 5 PTD edge extraction.
# GPU-side CSR adjacency will be built in Phase 5 from this dict if needed.
"""
struct TriMesh{T<:AbstractFloat, A<:AbstractMatrix{T}}
    vertices   :: A
    faces      :: Matrix{Int}
    normals    :: A
    areas      :: Vector{T}
    centroids  :: A
    edge_faces :: Dict{Tuple{Int,Int},Vector{Int}}
end

"""Number of vertices in `mesh`."""
nvertices(mesh::TriMesh) = size(mesh.vertices, 2)

"""Number of faces in `mesh`."""
nfaces(mesh::TriMesh) = size(mesh.faces, 2)

"""Signed surface area (total area of all faces) in `mesh`."""
total_area(mesh::TriMesh) = sum(mesh.areas)

# ── Internal helpers ──────────────────────────────────────────────────────────

# Per-face normals (unit), areas, centroids from vertex positions + face indices.
function _build_face_data(vertices::AbstractMatrix{T}, faces::Matrix{Int}) where {T<:AbstractFloat}
    Nf        = size(faces, 2)
    normals   = Matrix{T}(undef, 3, Nf)
    areas     = Vector{T}(undef, Nf)
    centroids = Matrix{T}(undef, 3, Nf)
    for i in 1:Nf
        i1, i2, i3 = faces[1, i], faces[2, i], faces[3, i]
        v1 = SVector{3,T}(vertices[1, i1], vertices[2, i1], vertices[3, i1])
        v2 = SVector{3,T}(vertices[1, i2], vertices[2, i2], vertices[3, i2])
        v3 = SVector{3,T}(vertices[1, i3], vertices[2, i3], vertices[3, i3])
        n  = cross(v2 - v1, v3 - v1)   # |n| = 2 * face area
        a  = norm(n) / 2
        areas[i]     = a
        # ponytail: degenerate face (a ≈ 0) gives zero normal; guard for robustness
        nhat = a > eps(T) ? n / (2a) : zero(SVector{3,T})
        normals[1, i] = nhat[1]; normals[2, i] = nhat[2]; normals[3, i] = nhat[3]
        c = (v1 + v2 + v3) / 3
        centroids[1, i] = c[1]; centroids[2, i] = c[2]; centroids[3, i] = c[3]
    end
    return normals, areas, centroids
end

# Signed volume via divergence theorem.
# Positive ↔ face normals point outward; negative ↔ inward; ≈0 for open/flat meshes.
function _signed_volume(vertices::AbstractMatrix{T}, faces::Matrix{Int}) where {T<:AbstractFloat}
    vol = zero(T)
    for i in 1:size(faces, 2)
        i1, i2, i3 = faces[1, i], faces[2, i], faces[3, i]
        v1 = SVector{3,T}(vertices[1, i1], vertices[2, i1], vertices[3, i1])
        v2 = SVector{3,T}(vertices[1, i2], vertices[2, i2], vertices[3, i2])
        v3 = SVector{3,T}(vertices[1, i3], vertices[2, i3], vertices[3, i3])
        vol += dot(v1, cross(v2, v3))
    end
    return vol / 6
end

# Undirected edge (sorted vertex-index pair) → list of incident face indices.
function _build_adjacency(faces::Matrix{Int})
    d = Dict{Tuple{Int,Int},Vector{Int}}()
    for i in 1:size(faces, 2)
        a, b, c = faces[1, i], faces[2, i], faces[3, i]
        for (u, v) in ((a, b), (b, c), (c, a))
            key = u < v ? (u, v) : (v, u)
            push!(get!(Vector{Int}, d, key), i)
        end
    end
    return d
end

# ── TriMesh constructor ───────────────────────────────────────────────────────

"""
    TriMesh(vertices, faces) -> TriMesh

Build a `TriMesh` from a 3×Nv vertex matrix and a 3×Nf face index matrix (1-based integers).

- Computes per-face normals, areas, and centroids.
- For closed meshes: if the signed volume (divergence-theorem) is negative, all face
  windings are flipped so normals point outward.  Open/flat meshes (signed volume ≈ 0)
  are left as-is — normals follow the initial winding.
- Builds edge-face adjacency for Phase 5 PTD extraction.
"""
function TriMesh(
    vertices :: AbstractMatrix{T},
    faces    :: AbstractMatrix{<:Integer},
) where {T<:AbstractFloat}
    F = convert(Matrix{Int}, faces)
    normals, areas, centroids = _build_face_data(vertices, F)
    if _signed_volume(vertices, F) < zero(T)
        # Swap rows 2 ↔ 3 reverses the cross-product → flips normal direction
        F[2, :], F[3, :] = copy(F[3, :]), copy(F[2, :])
        normals .= .-normals
    end
    edge_faces = _build_adjacency(F)
    return TriMesh{T, typeof(vertices)}(
        vertices, F, normals, areas, centroids, edge_faces,
    )
end

# ── I/O ──────────────────────────────────────────────────────────────────────

"""
    load_mesh(path, T=Float64) -> TriMesh

Load an STL or OBJ mesh via FileIO/MeshIO and return a `TriMesh{T}`.
Duplicate vertices (common in STL, which stores 3 independent verts per triangle)
are merged by exact position match before building adjacency.

# ponytail: exact-float deduplication (identical bit patterns).
# Upgrade to tolerance-based (e.g. k-d tree snap) if meshes arrive with
# floating-point jitter on nominally-shared vertices.
"""
function load_mesh(
    path :: AbstractString,
    ::Type{T} = Float64,
) where {T<:AbstractFloat}
    gb      = FileIO.load(path)
    raw_pts = GeometryBasics.coordinates(gb)   # Vector{Point{3,Tfloat}}
    raw_fcs = GeometryBasics.faces(gb)          # Vector{NgonFace{3,OffsetInteger{-1,UInt32}}}

    # Deduplicate vertices by exact position: build raw-index → dedup-index map.
    vertex_map   = Dict{NTuple{3,T},Int}()
    raw_to_dedup = Vector{Int}(undef, length(raw_pts))
    verts_list   = NTuple{3,T}[]
    for (idx, p) in enumerate(raw_pts)
        key      = (T(p[1]), T(p[2]), T(p[3]))
        dedup_i  = get!(vertex_map, key) do
            push!(verts_list, key)
            length(verts_list)
        end
        raw_to_dedup[idx] = dedup_i
    end

    Nv    = length(verts_list)
    Nf    = length(raw_fcs)
    verts = Matrix{T}(undef, 3, Nv)
    for (i, (x, y, z)) in enumerate(verts_list)
        verts[1, i] = x; verts[2, i] = y; verts[3, i] = z
    end
    faces_mat = Matrix{Int}(undef, 3, Nf)
    for (i, f) in enumerate(raw_fcs)
        for k in 1:3
            # convert(Int, OffsetInteger{-1, UInt32}) returns the 1-based index
            faces_mat[k, i] = raw_to_dedup[convert(Int, f[k])]
        end
    end

    return TriMesh(verts, faces_mat)
end

# ── Primitive builders ────────────────────────────────────────────────────────

"""
    unit_cube(; side=1.0, T=Float64) -> TriMesh

Axis-aligned cube `[0, side]³` as 12 triangles (2 per face, 6 faces).
Total surface area = 6·side².  Outward normals are exact axis directions (±x, ±y, ±z).
"""
function unit_cube(; side::Real = 1.0, T::Type{<:AbstractFloat} = Float64)
    s = T(side)
    # Vertices: corners of [0,s]^3  (3×8, one column per vertex)
    verts = Matrix{T}(undef, 3, 8)
    verts[:, 1] = [zero(T), zero(T), zero(T)]   # 1: (0,0,0)
    verts[:, 2] = [s,        zero(T), zero(T)]   # 2: (s,0,0)
    verts[:, 3] = [s,        s,        zero(T)]  # 3: (s,s,0)
    verts[:, 4] = [zero(T), s,        zero(T)]   # 4: (0,s,0)
    verts[:, 5] = [zero(T), zero(T), s       ]   # 5: (0,0,s)
    verts[:, 6] = [s,        zero(T), s       ]  # 6: (s,0,s)
    verts[:, 7] = [s,        s,        s      ]  # 7: (s,s,s)
    verts[:, 8] = [zero(T), s,        s       ]  # 8: (0,s,s)
    # 12 faces with CCW winding from outside → outward normals.
    # Each row = one vertex index of the face (column); 3 columns per face.
    # Verified by cross product: e.g. face (2,3,7): (v3-v2)×(v7-v2) ∝ +x ✓
    # Columns:   bot  bot  top  top  frt  frt  bck  bck  lft  lft  rgt  rgt
    faces = Int[1    1    5    5    1    1    4    4    1    1    2    2  ;
                4    3    6    7    2    6    8    7    5    8    3    7  ;
                3    2    7    8    6    5    7    3    8    4    7    6  ]
    return TriMesh(verts, faces)
end

# ── icosphere ─────────────────────────────────────────────────────────────────

# Add midpoint of edge (a, b) to verts_list, normalised to radius R (sphere projection).
# Shared between adjacent faces via midpoint_cache.
function _ico_midpoint!(
    verts_list     :: Vector{SVector{3,T}},
    midpoint_cache :: Dict{Tuple{Int,Int},Int},
    a :: Int, b :: Int, R :: T,
) where {T}
    key = a < b ? (a, b) : (b, a)
    return get!(midpoint_cache, key) do
        va = verts_list[a]
        vb = verts_list[b]
        mp = (va + vb) / 2
        push!(verts_list, mp * (R / norm(mp)))  # project to sphere
        length(verts_list)
    end
end

"""
    icosphere(n_subdiv=0; R=1.0, T=Float64) -> TriMesh

Regular icosphere: a regular icosahedron with `n_subdiv` rounds of midpoint
subdivision, each new vertex projected to the sphere of radius `R`.

| n_subdiv | faces | vertices |
|----------|-------|---------|
|    0     |    20 |      12  |
|    1     |    80 |      42  |
|    2     |   320 |     162  |
|    3     |  1280 |     642  |

Total surface area converges to 4π·R² from below as n_subdiv → ∞.

# ponytail: midpoint subdivision (not Loop or Catmull-Clark); no tangential
# smoothing is needed for Physical Optics because only normals/areas enter the
# PO integral, not mesh curvature.  Upgrade if the phase integrand needs
# smooth second derivatives (e.g. PTD correction in Phase 5).
"""
function icosphere(
    n_subdiv :: Int = 0;
    R        :: Real = 1.0,
    T        :: Type{<:AbstractFloat} = Float64,
)
    n_subdiv ≥ 0 || throw(ArgumentError("n_subdiv must be ≥ 0"))
    Rf = T(R)

    # Icosahedron base vertices (GLUT/standard): 12 vertices on unit sphere.
    # X = 1/√(1+φ²),  Z = φ/√(1+φ²)  where φ = (1+√5)/2.
    X = T(0.525731112119133606)
    Z = T(0.850650808352039932)
    raw_verts = SVector{3,T}[
        (-X,  0,  Z), ( X,  0,  Z), (-X,  0, -Z), ( X,  0, -Z),
        ( 0,  Z,  X), ( 0,  Z, -X), ( 0, -Z,  X), ( 0, -Z, -X),
        ( Z,  X,  0), (-Z,  X,  0), ( Z, -X,  0), (-Z, -X,  0),
    ]
    verts_list = [SVector{3,T}(v) * (Rf / norm(v)) for v in raw_verts]

    # 20 faces (GLUT 0-indexed → 1-indexed by +1).
    # Winding is consistently inward (signed volume < 0).
    # The TriMesh constructor's orient_outward! step flips them to outward.
    faces_list = NTuple{3,Int}[
        (1,5,2), (1,10,5), (10,6,5), (5,6,9), (5,9,2),
        (9,11,2),(9,4,11), (6,4,9),  (6,3,4), (3,8,4),
        (8,11,4),(8,7,11), (8,12,7), (12,1,7),(1,2,7),
        (7,2,11),(10,1,12),(10,12,3),(10,3,6),(8,3,12),
    ]

    for _ in 1:n_subdiv
        midpoint_cache = Dict{Tuple{Int,Int},Int}()
        new_faces = NTuple{3,Int}[]
        sizehint!(new_faces, 4 * length(faces_list))
        for (v1, v2, v3) in faces_list
            m12 = _ico_midpoint!(verts_list, midpoint_cache, v1, v2, Rf)
            m23 = _ico_midpoint!(verts_list, midpoint_cache, v2, v3, Rf)
            m13 = _ico_midpoint!(verts_list, midpoint_cache, v1, v3, Rf)
            push!(new_faces, (v1, m12, m13))
            push!(new_faces, (m12, v2, m23))
            push!(new_faces, (m13, m23, v3))
            push!(new_faces, (m12, m23, m13))
        end
        faces_list = new_faces
    end

    Nv    = length(verts_list)
    Nf    = length(faces_list)
    verts = Matrix{T}(undef, 3, Nv)
    for (i, v) in enumerate(verts_list)
        verts[1, i] = v[1]; verts[2, i] = v[2]; verts[3, i] = v[3]
    end
    faces_mat = Matrix{Int}(undef, 3, Nf)
    for (i, (a, b, c)) in enumerate(faces_list)
        faces_mat[1, i] = a; faces_mat[2, i] = b; faces_mat[3, i] = c
    end

    return TriMesh(verts, faces_mat)
end

# ── Flat plate ────────────────────────────────────────────────────────────────

"""
    flat_plate(Lx=1.0, Ly=1.0, nx=1, ny=1; T=Float64) -> TriMesh

Rectangular plate `[0,Lx]×[0,Ly]` in the z=0 plane, tessellated into 2·nx·ny
triangles.  Face normals point in the +z direction.

Not a closed mesh: `orient_outward!` is skipped (signed volume = 0).
"""
function flat_plate(
    Lx :: Real = 1.0, Ly :: Real = 1.0,
    nx :: Int  = 1,   ny :: Int  = 1;
    T  :: Type{<:AbstractFloat} = Float64,
)
    xs  = range(T(0), T(Lx); length = nx + 1)
    ys  = range(T(0), T(Ly); length = ny + 1)
    Nv  = (nx + 1) * (ny + 1)
    Nf  = 2 * nx * ny
    verts = Matrix{T}(undef, 3, Nv)
    faces = Matrix{Int}(undef, 3, Nf)

    # idx: column-major vertex index (i=x-step, j=y-step, both 1-based)
    idx(i, j) = (i - 1) * (ny + 1) + j

    k = 0
    for i in 1:nx+1, j in 1:ny+1
        k += 1
        verts[1, k] = xs[i]; verts[2, k] = ys[j]; verts[3, k] = zero(T)
    end

    k = 0
    for i in 1:nx, j in 1:ny
        v1 = idx(i, j);   v2 = idx(i+1, j)
        v3 = idx(i+1, j+1); v4 = idx(i, j+1)
        # CCW from +z → normal = +z
        k += 1; faces[:, k] = [v1, v2, v3]
        k += 1; faces[:, k] = [v1, v3, v4]
    end

    return TriMesh(verts, faces)
end

# ── Refinement ────────────────────────────────────────────────────────────────

"""
    refine(mesh::TriMesh) -> TriMesh

Uniform midpoint subdivision: each triangle → 4 by inserting edge midpoints.
Midpoints are shared between adjacent faces (no duplicate vertices on shared edges).
Face count is multiplied by exactly 4.

# ponytail: flat midpoint subdivision (not Loop).  Loop adds tangential vertex
# relaxation for smooth curvature, useful for visualisation but not required for
# Physical Optics where only areas and normals enter the scattering integral.
# Upgrade when the PTD correction (Phase 5) needs reliable curvature estimates.
"""
function refine(mesh::TriMesh{T}) where {T}
    verts   = mesh.vertices
    F       = mesh.faces
    Nv_old  = size(verts, 2)
    Nf_old  = size(F, 2)

    edge_mid  = Dict{Tuple{Int,Int},Int}()
    new_verts = SVector{3,T}[]
    sizehint!(new_verts, div(Nf_old * 3, 2))   # ≈ 1.5 × Nf unique edges

    midpoint_idx = let verts = verts, edge_mid = edge_mid, new_verts = new_verts, Nv_old = Nv_old
        function (a::Int, b::Int)
            key = a < b ? (a, b) : (b, a)
            get!(edge_mid, key) do
                va = SVector{3,T}(verts[1, a], verts[2, a], verts[3, a])
                vb = SVector{3,T}(verts[1, b], verts[2, b], verts[3, b])
                push!(new_verts, (va + vb) / 2)
                Nv_old + length(new_verts)
            end
        end
    end

    new_faces = Matrix{Int}(undef, 3, 4 * Nf_old)
    for i in 1:Nf_old
        v1, v2, v3 = F[1, i], F[2, i], F[3, i]
        m12 = midpoint_idx(v1, v2)
        m23 = midpoint_idx(v2, v3)
        m13 = midpoint_idx(v1, v3)
        new_faces[:, 4i-3] = [v1,  m12, m13]
        new_faces[:, 4i-2] = [m12, v2,  m23]
        new_faces[:, 4i-1] = [m13, m23, v3 ]
        new_faces[:, 4i  ] = [m12, m23, m13]
    end

    Nv_new        = Nv_old + length(new_verts)
    new_verts_mat = Matrix{T}(undef, 3, Nv_new)
    new_verts_mat[:, 1:Nv_old] .= verts
    for (i, v) in enumerate(new_verts)
        new_verts_mat[1, Nv_old+i] = v[1]
        new_verts_mat[2, Nv_old+i] = v[2]
        new_verts_mat[3, Nv_old+i] = v[3]
    end

    return TriMesh(new_verts_mat, new_faces)
end

"""
    refine(mesh::TriMesh, face_mask::AbstractVector{Bool}) -> TriMesh

Adaptive refinement: subdivide only the faces where `face_mask[i] == true`.
Unmarked faces are kept unchanged.

⚠ T-junctions: edges shared between a marked and an unmarked face will have a
midpoint vertex on the marked side but NOT on the unmarked side, creating a crack.
For crack-free results, propagate the mask to adjacent faces before calling or
run a second `refine` pass on the gap faces afterward.

# ponytail: T-junction allowed for Phase 3 GPU PO (the integrand is per-face,
# cracks do not affect the integral).  If Phase 5 PTD edge extraction requires
# watertight topology, add a gap-fill pass here.
"""
function refine(mesh::TriMesh{T}, face_mask::AbstractVector{Bool}) where {T}
    length(face_mask) == nfaces(mesh) ||
        throw(ArgumentError("face_mask length $(length(face_mask)) ≠ nfaces $(nfaces(mesh))"))

    verts     = mesh.vertices
    F         = mesh.faces
    Nv_old    = size(verts, 2)
    Nf_old    = size(F, 2)
    n_marked  = count(face_mask)

    edge_mid  = Dict{Tuple{Int,Int},Int}()
    new_verts = SVector{3,T}[]
    sizehint!(new_verts, div(n_marked * 3, 2))

    midpoint_idx = let verts = verts, edge_mid = edge_mid, new_verts = new_verts, Nv_old = Nv_old
        function (a::Int, b::Int)
            key = a < b ? (a, b) : (b, a)
            get!(edge_mid, key) do
                va = SVector{3,T}(verts[1, a], verts[2, a], verts[3, a])
                vb = SVector{3,T}(verts[1, b], verts[2, b], verts[3, b])
                push!(new_verts, (va + vb) / 2)
                Nv_old + length(new_verts)
            end
        end
    end

    Nf_new    = Nf_old + 3 * n_marked   # each marked face → +3 extra faces
    new_faces = Matrix{Int}(undef, 3, Nf_new)
    out       = 0
    for i in 1:Nf_old
        v1, v2, v3 = F[1, i], F[2, i], F[3, i]
        if face_mask[i]
            m12 = midpoint_idx(v1, v2)
            m23 = midpoint_idx(v2, v3)
            m13 = midpoint_idx(v1, v3)
            new_faces[:, out+1] = [v1,  m12, m13]
            new_faces[:, out+2] = [m12, v2,  m23]
            new_faces[:, out+3] = [m13, m23, v3 ]
            new_faces[:, out+4] = [m12, m23, m13]
            out += 4
        else
            new_faces[:, out+1] = [v1, v2, v3]
            out += 1
        end
    end

    Nv_new        = Nv_old + length(new_verts)
    new_verts_mat = Matrix{T}(undef, 3, Nv_new)
    new_verts_mat[:, 1:Nv_old] .= verts
    for (i, v) in enumerate(new_verts)
        new_verts_mat[1, Nv_old+i] = v[1]
        new_verts_mat[2, Nv_old+i] = v[2]
        new_verts_mat[3, Nv_old+i] = v[3]
    end

    return TriMesh(new_verts_mat, new_faces)
end

# ── Plotting stub ─────────────────────────────────────────────────────────────

"""
    plot_mesh(mesh::TriMesh; show_normals=false) -> Figure

Render the triangle mesh as a shaded 3-D surface in a new Makie Figure.
Requires a Makie backend to be loaded first, e.g. `using CairoMakie`.

With `show_normals=true`, arrows are drawn from each face centroid in the unit
outward-normal direction, scaled by √(mean face area) for visibility.
"""
function plot_mesh(mesh::TriMesh; kwargs...)
    ext = Base.get_extension(LowObservables, :LowObservablesMakieExt)
    ext === nothing &&
        error("Load a Makie backend first, e.g. `using CairoMakie`")
    return ext.plot_mesh(mesh; kwargs...)
end
