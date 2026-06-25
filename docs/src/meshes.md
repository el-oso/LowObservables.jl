# Surface Triangle Meshes

This page is for readers who are new to 3-D mesh representations.
If you already know what a surface mesh is, jump to the [API Reference](api.md).

---

## Why triangles?

A real radar target — a jet, a missile body, a ship — has a complex curved surface.
A computer cannot represent a smooth curve exactly; it approximates it with
a **polygon mesh**: a collection of flat polygons that tile the surface.

Triangles are the simplest polygon.  Every triangle is guaranteed to be:

- **Planar** — its three vertices always lie in a single plane (a property that
  quadrilaterals can lose if the fourth vertex drifts).
- **Convex** — no concave angles.
- **Well-defined** — no ambiguity about which way the face "points".

For Physical Optics (Phase 3), each triangle becomes one panel in a quadrature
sum.  The finer the mesh, the better the sum approximates the true surface integral.

---

## What a TriMesh stores

`LowObservables.TriMesh{T, A}` is a **struct-of-arrays** mesh.

### Vertices

An `A`-matrix of size `3 × Nv` where column `k` holds the `(x, y, z)` position
of vertex `k` in metres.  `A` defaults to `Matrix{T}` on CPU; Phase 3 will
accept `CuMatrix{T}` for GPU kernels.

### Faces

A `3 × Nf` integer matrix.  Column `k` holds the three 1-based vertex indices
of the k-th triangle.  Indices are always consistent: the vertices are listed
in **counter-clockwise order when viewed from the outside**, so the face's
outward normal follows the right-hand rule.

### Per-face geometry

These are what the PO kernel reads for each face:

| Field        | Size     | Meaning                                  |
|--------------|----------|------------------------------------------|
| `normals`    | 3 × Nf   | Unit outward normal vector               |
| `areas`      | Nf       | Face area [m²]                           |
| `centroids`  | 3 × Nf   | Face centroid (arithmetic mean of verts) |

---

## Outward normals

For a **closed** mesh (like a sphere or a cube), every face has an inside and an
outside.  The Physical Optics integral needs the outward normal — the one pointing
away from the interior.

LowObservables checks the orientation automatically using the **divergence theorem**:
the signed volume of a closed mesh is positive when normals point outward.
If the computed signed volume is negative (inward winding), all face windings
are flipped.

---

## Refinement: why it matters

The Physical Optics integral approximates the surface field as constant across
each face.  The error is proportional to the face area:

- Coarse mesh (few, large triangles) → fast but inaccurate.
- Fine mesh (many, small triangles) → slow but converges to the exact integral.

`refine(mesh)` splits every triangle into **four smaller triangles** by inserting
edge midpoints.  Each call multiplies the face count by 4:

| Calls to `refine` | Faces (icosphere L0 = 20) |
|-------------------|---------------------------|
| 0                 | 20                        |
| 1                 | 80                        |
| 2                 | 320                       |
| 3                 | 1 280                     |
| 4                 | 5 120                     |

For a **flat** mesh, midpoint subdivision preserves the total area exactly.
For a sphere approximated by an icosahedron, the icosphere primitive projects
new vertices back onto the sphere at each subdivision level — so each level
gets closer to the true surface area 4πR².

---

## Interactive: icosphere refinement

The widget below shows how successive refinements of an icosahedron approximate
the unit sphere.  The left panel shows the shaded surface; the right panel shows
per-face outward normals as red arrows.

**Drag the slider** to change the refinement level (0 = icosahedron, 3 = 1280 faces)
and watch the surface become smoother.  The widget runs entirely in your browser.

```@setup meshes
using WGLMakie, Bonito
WGLMakie.activate!()
Makie.inline!(true)
Page(exportable = true, offline = true)
```

```@example meshes
using LowObservables, WGLMakie, Bonito

# Precompute icospheres at refinement levels 0..3 and convert to GeometryBasics
# meshes for Makie's mesh! recipe.
# ponytail: precompute all levels; record_states captures each slider state
# so the widget works offline.  Live re-subdivision per tick is ~seconds per
# frame on a 1280-face mesh.
n_levels = 4
lo_meshes = [icosphere(n) for n in 0:n_levels-1]

function _to_gb(m)
    pts  = [WGLMakie.Point3f(Float32(m.vertices[1,i]),
                              Float32(m.vertices[2,i]),
                              Float32(m.vertices[3,i])) for i in 1:nvertices(m)]
    tris = [WGLMakie.GLTriangleFace(m.faces[1,j], m.faces[2,j], m.faces[3,j])
            for j in 1:nfaces(m)]
    Makie.GeometryBasics.Mesh(pts, tris)
end

gb_meshes = [_to_gb(m) for m in lo_meshes]

App() do session::Session
    sl = Bonito.Slider(1:n_levels)

    fig = Figure(size = (680, 760))
    ax1 = Axis3(fig[1, 1]; title = "Icosphere — shaded surface", aspect = :equal)
    ax2 = Axis3(fig[2, 1]; title = "Per-face normals", aspect = :equal)

    # Observable mesh driven by slider
    obs_mesh = map(sl.value) do lv; gb_meshes[lv] end
    mesh!(ax1, obs_mesh; color = :steelblue, shading = FastShading)

    # Observable normal arrows (scale by √mean_area for visibility)
    obs_arrow_pts = map(sl.value) do lv
        m = lo_meshes[lv]
        [WGLMakie.Point3f(Float32(m.centroids[1,i]),
                           Float32(m.centroids[2,i]),
                           Float32(m.centroids[3,i])) for i in 1:nfaces(m)]
    end
    obs_arrow_dirs = map(sl.value) do lv
        m = lo_meshes[lv]
        sc = Float32(sqrt(sum(m.areas) / nfaces(m)) * 0.5)
        [WGLMakie.Vec3f(Float32(m.normals[1,i]),
                         Float32(m.normals[2,i]),
                         Float32(m.normals[3,i])) * sc for i in 1:nfaces(m)]
    end
    mesh!(ax2, obs_mesh; color = (:steelblue, 0.3), shading = FastShading)
    arrows3d!(ax2, obs_arrow_pts, obs_arrow_dirs; color = :tomato, shaftradius = 0.01)

    # Stats label
    label = map(sl.value) do lv
        m = lo_meshes[lv]
        a = total_area(m)
        "Level $(lv-1): $(nfaces(m)) faces, $(nvertices(m)) verts  |  area = $(round(a, digits=4)) m²  (4πR² = $(round(4π, digits=4)))"
    end

    ui = DOM.div("Refinement level: ", sl, DOM.br(), label)
    return Bonito.record_states(session, DOM.div(ui, fig))
end
```

---

## Primitive shapes

```julia
using LowObservables

# Axis-aligned unit cube (12 triangles, area = 6)
cube = unit_cube()

# Icosphere with 2 rounds of subdivision (320 faces)
sphere = icosphere(2; R = 1.0)

# Flat rectangular plate, 4×4 quads (32 triangles)
plate = flat_plate(2.0, 1.0, 4, 4)

# Load an STL or OBJ file
# mesh = load_mesh("my_target.stl")
```

## Refine and visualise

```julia
using LowObservables, CairoMakie

s = icosphere(1)               # 80 faces
s2 = refine(s)                 # 320 faces — flat midpoint subdivision

fig = plot_mesh(s2; show_normals = true)   # requires a Makie backend
```

`refine(mesh, mask)` subdivides only the faces where `mask[i] == true`,
enabling adaptive refinement where the scattering integrand varies rapidly.
