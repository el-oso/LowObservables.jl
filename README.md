# LowObservables.jl

[![CI](https://github.com/el-oso/LowObservables.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/el-oso/LowObservables.jl/actions/workflows/CI.yml)
[![Coverage](https://coveralls.io/repos/github/el-oso/LowObservables.jl/badge.svg?branch=master)](https://coveralls.io/github/el-oso/LowObservables.jl?branch=master)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Radar cross-section (RCS) simulation for Julia.

**Phase 0** ships the analytically-exact Mie series for a perfectly-conducting sphere —
the standard validation target for every RCS code.

## Features

- **Exact PEC-sphere Mie series** — monostatic (backscatter) RCS σ [m²] via the
  Bohren & Huffman formulation with the Wiscombe (1980) truncation criterion
- **Generic float type** — `mie_rcs_sphere(Float32(a), Float32(λ))` returns `Float32`;
  ready for GPU arrays and automatic differentiation in later phases
- **TypeContracts interface** — `AbstractRCSModel` contract; `@verify`-checked at
  load time; extensible for GTD/PO models in later phases
- **Validated numerics** — tests cover Rayleigh ∝ (ka)⁴ limit, optical limit πa² (within 0.03%), and first-resonance peak σ/πa² ≈ 3.65

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/el-oso/LowObservables.jl")
```

## Quick start

```julia
using LowObservables

a = 0.10        # sphere radius [m]
λ = 0.03        # wavelength [m]  (≈ 10 GHz)
σ = mie_rcs_sphere(a, λ)              # → ~0.030 m²
σ_dBsm = 10log10(σ)                   # → ~ −15 dBsm

# sweep over wavelengths
λs  = range(0.001, 0.3, length=500)
σs  = rcs_vs_size(a, λs)
σ_opt = optical_limit_rcs(a)          # πa² = 0.0314 m²

# model-type API
m = MieSphere(a)
rcs(m, λ)                              # same result
```

## Validation

| `ka` | σ/(πa²) | notes |
|------|---------|-------|
| 0.01 | 9.0×10⁻⁸ | Rayleigh: 9(ka)⁴ match to <0.002% |
| 1.03 | 3.655 | first resonance peak (literature ≈ 3.6–4.0) |
| [30,31] mean | 0.9997 | optical limit (0.03% error) |

## License

MIT — see [LICENSE](LICENSE).
