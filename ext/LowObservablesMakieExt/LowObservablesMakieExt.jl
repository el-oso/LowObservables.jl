"""
Makie rendering for LowObservables triangle meshes.

Loaded automatically when any Makie backend (CairoMakie, GLMakie, WGLMakie)
is imported alongside LowObservables.  No direct dep on GeometryBasics —
use `Makie.GeometryBasics` for GB types since GB is a transitive dep only.
"""
module LowObservablesMakieExt

using Makie
using LowObservables
using LowObservables: TriMesh, nvertices, nfaces

export plot_mesh, plot_mesh!

# Convert a TriMesh to a GeometryBasics.Mesh for Makie's mesh! recipe.
function _to_gb_mesh(mesh::TriMesh{T}) where {T}
    Nv   = nvertices(mesh)
    Nf   = nfaces(mesh)
    pts  = [Makie.GeometryBasics.Point3f(
                Float32(mesh.vertices[1, i]),
                Float32(mesh.vertices[2, i]),
                Float32(mesh.vertices[3, i]),
            ) for i in 1:Nv]
    tris = [Makie.GeometryBasics.GLTriangleFace(
                mesh.faces[1, j],
                mesh.faces[2, j],
                mesh.faces[3, j],
            ) for j in 1:Nf]
    return Makie.GeometryBasics.Mesh(pts, tris)
end

"""
    plot_mesh!(ax::Axis3, mesh::TriMesh; show_normals=false, kw...) -> Axis3

Render `mesh` into an existing `Axis3` as a shaded surface.
With `show_normals=true`, arrows are drawn from each face centroid in the outward
normal direction, scaled by √(mean face area) for visibility.

Returns `ax` for chaining.
"""
function plot_mesh!(
    ax          :: Makie.Axis3,
    mesh        :: TriMesh{T};
    show_normals :: Bool = false,
    color        = :steelblue,
    kwargs...
) where {T}
    gb = _to_gb_mesh(mesh)
    Makie.mesh!(ax, gb; color, kwargs...)

    if show_normals
        Nf        = nfaces(mesh)
        mean_area = sum(mesh.areas) / Nf
        scale     = Float32(sqrt(mean_area) * 0.4)
        pts   = [Makie.Point3f(
                     Float32(mesh.centroids[1, i]),
                     Float32(mesh.centroids[2, i]),
                     Float32(mesh.centroids[3, i]),
                 ) for i in 1:Nf]
        dirs  = [Makie.Vec3f(
                     Float32(mesh.normals[1, i]),
                     Float32(mesh.normals[2, i]),
                     Float32(mesh.normals[3, i]),
                 ) * scale for i in 1:Nf]
        # ponytail: arrows3d! (Makie 0.24+); arrows! is deprecated in 0.24
        Makie.arrows3d!(ax, pts, dirs; color = :tomato, shaftradius = 0.02 * scale)
    end
    return ax
end

"""
    plot_mesh(mesh::TriMesh; show_normals=false, kw...) -> Figure

Render `mesh` in a new Makie Figure.
With `show_normals=true`, outward normal arrows are drawn from face centroids.

```julia
using LowObservables, CairoMakie
mesh = icosphere(2)
fig  = plot_mesh(mesh; show_normals=true)
```
"""
function plot_mesh(mesh::TriMesh; kwargs...)
    fig = Makie.Figure()
    ax  = Makie.Axis3(fig[1, 1]; aspect = :equal)
    plot_mesh!(ax, mesh; kwargs...)
    return fig
end

end # module
