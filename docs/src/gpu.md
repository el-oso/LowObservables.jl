# GPU acceleration (CUDA)

The Physical Optics solver is written with
[KernelAbstractions](https://github.com/JuliaGPU/KernelAbstractions.jl): the kernel selects
its hardware backend from the arrays it is given (`get_backend(normals)`). So running on a GPU
needs **no kernel changes** — only the mesh's per-face arrays placed on the device. That is what
[`to_gpu`](@ref) does (it requires `CUDA.jl` to be loaded; the conversion lives in a package
extension, so CUDA is a *weak* dependency and CPU-only installs never pull it in).

```julia
using LowObservables, CUDA

mesh = icosphere(5; R = 1.0)          # 20 480 facets
k    = 2π / 0.1
ki   = [0.0, 0.0, 1.0]; ei = [1.0, 0.0, 0.0]

σ_cpu = po_rcs_monostatic(mesh,         ki, ei; k = k)   # CPU
σ_gpu = po_rcs_monostatic(to_gpu(mesh), ki, ei; k = k)   # same call, runs on the GPU
```

## Correctness

The GPU result matches the CPU result to machine precision — the same kernel, just a different
backend. On an RTX 4070:

| mesh | facets | CPU σ | GPU σ |
|------|--------|-------|-------|
| icosphere(3) | 1 280 | 98.6236 | 98.6236 |
| icosphere(4) | 5 120 | 12.3827 | 12.3827 |
| icosphere(5) | 20 480 | 3.06416 | 3.06416 |

An aspect sweep agrees to `max reldiff ≈ 3×10⁻¹⁴`.

## Performance

A 180-direction monostatic sweep over the 20 480-facet mesh (RTX 4070):

| | time |
|---|---|
| CPU | 53.2 ms |
| GPU | 30.3 ms |
| **speedup** | **1.8×** |

This is a modest, honest win at this scale: 20 k facets do not saturate the GPU, and the sweep
currently loops directions on the host — one kernel launch plus a reduction per direction — so
launch and reduction overhead dominate. The speedup grows with mesh size and would improve
substantially with a kernel batched over all directions at once (a natural next step).

!!! note "Tested on hardware"
    The CUDA path is validated on an NVIDIA RTX 4070 (the [`cuda_test.jl`](https://github.com/el-oso/LowObservables.jl/blob/master/test/cuda_test.jl)
    test runs there and **skips** automatically on CPU-only machines and CI). Benchmark datapoints
    are saved under `bench/results/`.
