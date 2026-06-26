"""
Shape optimization for radar cross-section reduction.

Minimizes the mean monostatic PO RCS over a threat sector by deforming mesh
vertices.  Gradient is computed via ForwardDiff (forward-mode AD) on a
plain-Julia differentiable PO objective function.  Optim.jl (LBFGS) performs
the descent.

Objective:  J(x) = sector_rcs(x) + λ_reg · ‖x − x₀‖²

The regulariser is a Tikhonov penalty on displacement from the initial vertex
positions x₀.  Without it, the optimizer collapses the mesh to a point (zero
area → zero RCS), which is physically meaningless.

Physical intuition (F-117-style shaping):
  A flat panel illuminated near-normally creates a specular flash toward the
  radar — the classic RCS spike.  Tilting the panel reflects the flash away
  from the threat direction.  The optimizer finds vertex displacements that
  tilt panels to deflect specular returns outside the threat sector.

# ponytail: ForwardDiff (forward-mode AD) costs O(3Nv) forward passes per
#   gradient evaluation.  Switch to Enzyme (reverse-mode, O(1)) for Nv >> 1000.
# ponytail: Float64 throughout; Optim's LBFGS line-search uses Float64 internally.
#   Float32 meshes are accepted but vertices are promoted to Float64 for the solve.
# ponytail: piecewise-constant illumination indicator (dot(ki,n)<0 → Bool).
#   Its derivative is zero a.e., which is the correct PO shape-derivative.
#   Facets that flip lit↔shadow between steps cause a J discontinuity but not a
#   gradient error — acceptable for gradient descent.

References:
  Knott et al., *Radar Cross Section*, 2nd ed. (2004) §9 (shaping for stealth)
  Ufimtsev (2014), §6.3 (PO shape derivatives)
"""

using ForwardDiff
using Optim

# ── Per-face differentiable PO contribution ───────────────────────────────────

# Compute the real/imag parts of the PO G-vector contribution from face j.
# x      — flat vertex vector (length 3Nv), column-major: vertex k at x[3k-2:3k]
# i1,i2,i3 — 1-based vertex indices for this face
# ki     — incident direction NTuple{3,Float64} (not differentiated)
# smki   — s − ki = −2ki (monostatic), NTuple{3,Float64}
# kice   — ki × ei, NTuple{3,Float64}
# k      — wavenumber Float64
#
# Dual-safe: the only branch is on the SIGN of dot(ki,n), which is piecewise-
# constant (returns Bool for both Float64 and ForwardDiff.Dual inputs, because
# ForwardDiff.Dual <: Real and isless compares value parts).
@inline function _face_G(x, i1, i2, i3, ki, smki, kice, k)
    # vertex coordinates (Dual when differentiating)
    v1x = x[3i1-2]; v1y = x[3i1-1]; v1z = x[3i1]
    v2x = x[3i2-2]; v2y = x[3i2-1]; v2z = x[3i2]
    v3x = x[3i3-2]; v3y = x[3i3-1]; v3z = x[3i3]

    # edge vectors
    e1x = v2x-v1x; e1y = v2y-v1y; e1z = v2z-v1z
    e2x = v3x-v1x; e2y = v3y-v1y; e2z = v3z-v1z

    # unnormalized cross product n_cross = e1 × e2  (‖n_cross‖ = 2·area)
    ncx = e1y*e2z - e1z*e2y
    ncy = e1z*e2x - e1x*e2z
    ncz = e1x*e2y - e1y*e2x
    nc_norm = sqrt(ncx^2 + ncy^2 + ncz^2)
    area    = nc_norm / 2

    # unit normal
    inv_nc = 1 / nc_norm
    nx = ncx * inv_nc; ny = ncy * inv_nc; nz = ncz * inv_nc

    # Illumination: branch only on sign — returns Bool for both Float64 and Dual.
    # ponytail: indicator is piecewise-constant; derivative is 0 a.e. (correct).
    if ki[1]*nx + ki[2]*ny + ki[3]*nz < 0
        # centroid
        cx = (v1x + v2x + v3x) / 3
        cy = (v1y + v2y + v3y) / 3
        cz = (v1z + v2z + v3z) / 3

        # n̂ × kice  (kice = ki × ei)
        bx = ny*kice[3] - nz*kice[2]
        by = nz*kice[1] - nx*kice[3]
        bz = nx*kice[2] - ny*kice[1]

        # phase = k · dot(s−ki, centroid)
        ph = k * (smki[1]*cx + smki[2]*cy + smki[3]*cz)
        sr = cos(ph); si = sin(ph)

        return (area*bx*sr, area*bx*si,
                area*by*sr, area*by*si,
                area*bz*sr, area*bz*si)
    else
        z = zero(area)   # zero(Dual{...}) is defined
        return (z, z, z, z, z, z)
    end
end

# ── Differentiable total surface area ─────────────────────────────────────────

# Sum of all face areas from the flat vertex vector x.  ForwardDiff-able.
# Used by the size-preservation constraint in optimize_rcs.
function _total_area(x, faces)
    A = zero(eltype(x))
    for j in axes(faces, 2)
        i1, i2, i3 = faces[1,j], faces[2,j], faces[3,j]
        e1x = x[3i2-2]-x[3i1-2]; e1y = x[3i2-1]-x[3i1-1]; e1z = x[3i2]-x[3i1]
        e2x = x[3i3-2]-x[3i1-2]; e2y = x[3i3-1]-x[3i1-1]; e2z = x[3i3]-x[3i1]
        ncx = e1y*e2z - e1z*e2y
        ncy = e1z*e2x - e1x*e2z
        ncz = e1x*e2y - e1y*e2x
        A += sqrt(ncx^2 + ncy^2 + ncz^2) / 2
    end
    return A
end

# ── Public: differentiable sector-averaged RCS ────────────────────────────────

"""
    sector_rcs(x, faces, dirs, ei; k) → σ̄ [m²]

Mean monostatic PO RCS over the threat sector `dirs` (3×N matrix of unit
incidence directions) for a mesh whose vertex positions are encoded in the flat
vector `x` (length 3·Nv, column-major reshape of the 3×Nv vertex matrix).

ForwardDiff-differentiable: all face geometry (centroid, area, unit normal) is
computed analytically from `x`.  The lit/shadow indicator `dot(ki,n)<0` is
piecewise-constant (its derivative is zero a.e.) and correctly handled by
ForwardDiff because `<` on Dual numbers compares value parts → returns Bool.

Arguments:
- `x`     — flat Float64 (or Dual) vertex vector, length 3·Nv
- `faces` — 3×Nf integer face-index matrix (fixed; not differentiated)
- `dirs`  — 3×N threat-sector incidence directions (normalised internally)
- `ei`    — polarisation unit vector (fixed; not differentiated)
- `k`     — wavenumber [rad/m]
"""
function sector_rcs(
    x     :: AbstractVector,
    faces :: Matrix{Int},
    dirs  :: AbstractMatrix,
    ei;
    k     :: Real,
)
    Nd = size(dirs, 2)
    Nf = size(faces, 2)
    kF = Float64(k)

    # Fixed (non-differentiated) unit vectors — plain Float64 tuples
    ei_n = let n = sqrt(ei[1]^2 + ei[2]^2 + ei[3]^2)
        (Float64(ei[1])/n, Float64(ei[2])/n, Float64(ei[3])/n)
    end

    σ_sum = zero(eltype(x))

    for d in 1:Nd
        # incident direction — plain Float64
        ki = let c = view(dirs, :, d), n = sqrt(c[1]^2 + c[2]^2 + c[3]^2)
            (Float64(c[1])/n, Float64(c[2])/n, Float64(c[3])/n)
        end
        s    = (-ki[1], -ki[2], -ki[3])   # monostatic: ŝ = −k̂_i
        kice = (ki[2]*ei_n[3] - ki[3]*ei_n[2],
                ki[3]*ei_n[1] - ki[1]*ei_n[3],
                ki[1]*ei_n[2] - ki[2]*ei_n[1])
        smki = (s[1]-ki[1], s[2]-ki[2], s[3]-ki[3])   # = −2ki

        # Accumulate G = Σ_lit A_i (n̂_i × kice) exp(i·phase_i)
        Gr1 = Gr2 = Gr3 = Gi1 = Gi2 = Gi3 = zero(eltype(x))
        for j in 1:Nf
            i1, i2, i3 = faces[1,j], faces[2,j], faces[3,j]
            gre1, gim1, gre2, gim2, gre3, gim3 =
                _face_G(x, i1, i2, i3, ki, smki, kice, kF)
            Gr1 += gre1;  Gi1 += gim1
            Gr2 += gre2;  Gi2 += gim2
            Gr3 += gre3;  Gi3 += gim3
        end

        # Gperp = G − dot(ŝ, G) · ŝ   (transverse component)
        dre = s[1]*Gr1 + s[2]*Gr2 + s[3]*Gr3
        dim = s[1]*Gi1 + s[2]*Gi2 + s[3]*Gi3

        # σ = (k²/π) · |Gperp|²  — written as sum of squares (avoids Complex{Dual})
        σ_sum += (kF^2 / π) * (
            (Gr1 - dre*s[1])^2 + (Gi1 - dim*s[1])^2 +
            (Gr2 - dre*s[2])^2 + (Gi2 - dim*s[2])^2 +
            (Gr3 - dre*s[3])^2 + (Gi3 - dim*s[3])^2
        )
    end

    return σ_sum / Nd
end

# ── Public: optimizer ─────────────────────────────────────────────────────────

"""
    optimize_rcs(mesh, dirs, ei; k, μ_area, λ_reg, maxiters) → (; mesh, history)

Deform `mesh` vertices to minimise the mean monostatic PO RCS over the threat
sector `dirs` (3×N direction matrix).

Objective:  `J(x) = sector_rcs(x) + μ_area·((A(x)−A₀)/A₀)² + λ_reg·‖x−x₀‖²`

The **size-preservation** term `μ_area·((A(x)−A₀)/A₀)²` holds the total surface
area near its initial value `A₀`.  Without it the optimizer takes the trivial
collapse solution — shrinking the mesh (σ ∝ A² for a plate) instead of
**reshaping** it.  With size held, the only way to lower σ is to *tilt/curve*
the surface so the specular flash is deflected out of the threat sector and the
facet phases interfere destructively.  This is the physical content of stealth
shaping.

The small Tikhonov term `λ_reg·‖x−x₀‖²` anchors the translation/rigid null modes
(which leave σ unchanged) so the solver stays well-conditioned.

**Symmetry breaking (`seed_eps`).**  A flat plate at normal incidence is a
*stationary point* of σ: by the reflection symmetry z → −z, σ is even in any
z-deformation, so the z-gradient is exactly zero at the flat start and gradient
descent never leaves the flat manifold (it would only find the now-penalised
shrink direction).  `seed_eps` adds a tiny deterministic perturbation to the
starting vertices to break this symmetry; the optimizer then rolls downhill into
the tilt/curve solution.  Set `seed_eps = 0` to disable (e.g. when the initial
shape is already asymmetric).

# ponytail: area constraint for OPEN meshes (the plate demo).  For a CLOSED mesh
#   swap in enclosed VOLUME preservation (divergence-theorem signed volume) so a
#   convex body cannot deflate.  Not wired up — add when optimizing closed shells.

Returns a named tuple:
- `mesh`    — optimised `TriMesh` (reconstructed from optimised vertices + original faces)
- `history` — `(; σ, vertices)`:
    - `σ`        :: `Vector{Float64}` — sector-mean RCS [m²] at each stored iteration
    - `vertices` :: `Vector{Matrix{Float64}}` — 3×Nv vertex snapshots (initial + each iter)

# ponytail: a perfectly symmetric start (flat plate at normal incidence) is a
#   BIFURCATION — folding up vs down deflects the specular equally, so two mirror
#   minima exist at identical depth.  Near that knife-edge the converged branch is
#   bistable to FP-level perturbations (a geometric degeneracy, not solver
#   nondeterminism).  Well-resolved problems have a deep basin and converge firmly;
#   if you need a guaranteed branch, start from a slightly non-planar mesh.
# ponytail: ForwardDiff forward-mode, O(3Nv) passes per gradient.
#   Use Enzyme reverse-mode for Nv ≫ 1000.
# ponytail: LBFGS (quasi-Newton); robust for smooth objectives.
#   Try GradientDescent if LBFGS diverges on highly non-convex shapes.
# ponytail: extended_trace stores per-iteration x copies; memory O(maxiters × 3Nv).
#   Disable for large meshes by passing record_history=false (not yet wired up).
"""
function optimize_rcs(
    mesh     :: TriMesh,
    dirs     :: AbstractMatrix,
    ei;
    k        :: Real    = 2π / 0.1,
    μ_area   :: Real    = 1e5,
    λ_reg    :: Real    = 1e-2,
    seed_eps :: Real    = 1e-3,
    maxiters :: Int     = 50,
)
    x0 = vec(Float64.(mesh.vertices))   # always Float64 for the optimizer
    F  = mesh.faces
    Nv = size(mesh.vertices, 2)
    kF = Float64(k)
    μF = Float64(μ_area)
    λF = Float64(λ_reg)
    A0 = _total_area(x0, F)

    # Deterministic symmetry-breaking seed (fixed hash → reproducible, no RNG dep).
    # Anchored by λ_reg back toward x0, so the seed only chooses a descent branch.
    # ponytail: GLSL-style frac(sin) hash; deterministic pseudo-noise in [-1,1].
    x_start = copy(x0)
    if seed_eps > 0
        εs = Float64(seed_eps)
        for i in eachindex(x_start)
            h = sin(i * 12.9898) * 43758.5453
            x_start[i] += εs * (2*(h - floor(h)) - 1)
        end
    end

    # Objective: RCS + strong size-preservation + weak displacement anchor
    function J(x)
        σ    = sector_rcs(x, F, dirs, ei; k = kF)
        area = _total_area(x, F)
        size_pen = μF * ((area - A0) / A0)^2
        anchor   = λF * sum((xi - x0i)^2 for (xi, x0i) in zip(x, x0))
        return σ + size_pen + anchor
    end

    # In-place gradient (ForwardDiff)
    function J_grad!(g, x)
        g .= ForwardDiff.gradient(J, x)
        return nothing
    end

    result = Optim.optimize(
        J, J_grad!, x_start, LBFGS(),
        Optim.Options(
            iterations     = maxiters,
            store_trace    = true,
            extended_trace = true,   # stores metadata["x"] per iteration
            show_trace     = false,
        ),
    )

    tr     = Optim.trace(result)
    x_opt  = Optim.minimizer(result)

    # Build history: initial state + one entry per recorded iteration
    all_x   = vcat([x0], [copy(s.metadata["x"]) for s in tr])
    hist_σ  = [sector_rcs(xi, F, dirs, ei; k = kF) for xi in all_x]
    hist_v  = [reshape(xi, 3, Nv) for xi in all_x]

    mesh_opt = TriMesh(reshape(x_opt, 3, Nv), F)
    return (; mesh = mesh_opt, history = (; σ = hist_σ, vertices = hist_v))
end
