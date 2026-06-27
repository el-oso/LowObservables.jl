# GPU (CUDA) backend tests. These run only where CUDA.jl is installed AND a device is
# present (e.g. the galen RTX 4070 box); they SKIP cleanly on CPU-only machines and CI so
# the suite stays green everywhere without forcing a heavy CUDA download.
#
# CUDA is a weak dependency, so we locate/require it by PkgId rather than `using CUDA`
# (which would need CUDA in test/Project.toml). On a GPU box: add CUDA to the test env
# (`Pkg.add("CUDA")`) and these run.

@testitem "CUDA backend: PO RCS matches CPU" begin
    using LowObservables

    cuda_id = Base.PkgId(Base.UUID("052768ef-5323-5732-b1bb-66c8b64840ba"), "CUDA")
    if Base.locate_package(cuda_id) === nothing
        @test_skip "CUDA.jl not installed in this environment"
    else
        CUDA = Base.require(cuda_id)           # loads CUDA → triggers LowObservablesCUDAExt
        if !CUDA.functional()
            @test_skip "no functional CUDA device"
        else
            k  = 2π / 0.1
            ei = [1.0, 0.0, 0.0]

            # monostatic: CPU vs GPU on the same mesh
            m   = icosphere(3; R = 1.0)
            ki  = [0.0, 0.0, 1.0]
            cpu = po_rcs_monostatic(m, ki, ei; k = k)
            gpu = po_rcs_monostatic(to_gpu(m), ki, ei; k = k)
            @test isapprox(cpu, gpu; rtol = 1e-4)

            # to_gpu actually moved the per-face arrays onto the device
            mg = to_gpu(m)
            @test !(mg.normals isa Array)
            @test occursin("CuArray", string(typeof(mg.areas)))

            # an aspect sweep agrees too
            θs   = range(0, π; length = 24)
            dirs = reduce(hcat, ([-sin(θ), 0.0, -cos(θ)] for θ in θs))
            cpus = po_rcs_sweep(m, dirs, ei; k = k)
            gpus = po_rcs_sweep(mg, dirs, ei; k = k)
            @test isapprox(cpus, gpus; rtol = 1e-4)
        end
    end
end
