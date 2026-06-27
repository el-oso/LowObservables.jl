module LowObservablesCUDAExt

# CUDA backend for the Physical Optics solver.
#
# The PO kernel (`src/physical_optics.jl`) is written with KernelAbstractions and selects
# its backend via `get_backend(normals)`. So all that is needed for GPU execution is to
# place the mesh's per-face arrays on the device — then the existing `po_rcs_*` functions
# run on CUDA unchanged. `to_gpu` does exactly that.

using LowObservables: LowObservables, TriMesh
using CUDA: CuArray

function LowObservables.to_gpu(m::TriMesh)
    return TriMesh(
        CuArray(m.vertices),   # 3×Nv
        m.faces,               # kept on CPU — the PO kernel never reads connectivity
        CuArray(m.normals),    # 3×Nf  (read by the kernel)
        CuArray(m.areas),      # Nf    (read by the kernel)
        CuArray(m.centroids),  # 3×Nf  (read by the kernel)
        m.edge_faces,          # CPU adjacency Dict (PTD preprocessing, not the PO kernel)
    )
end

end
