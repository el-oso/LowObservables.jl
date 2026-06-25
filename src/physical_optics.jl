"""
Physical Optics (PO) monostatic/bistatic RCS of a PEC triangulated surface.

Convention: e^{+jωt}.  Plane wave propagates in direction `ki` (unit vector),
scattered toward `s` (unit vector), incident polarisation `ei` (unit, ⊥ ki).
Wavenumber k = 2π/λ.  Face i is lit iff dot(ki, n_i) < 0.

Surface-current model: J = 2 n̂ × H on the lit side, J = 0 on the shadow side.
This is the Physical Optics (Kirchhoff) approximation — exact at high frequency
for smooth surfaces, breaks near shadow boundaries and edges.

Per-facet centroid quadrature:

    G = Σ_lit  A_i (n_i × (ki × ei)) exp(im·k·dot(s−ki, c_i))
    Gperp = G − dot(s, G)·s        (transverse component; s real, no conjugation)
    σ = (k²/π) |Gperp|²_vec        (sum of abs² of complex components)

Monostatic: s = −ki  →  phase factor exp(−im·2k·dot(ki, c_i)).

# ponytail: centroid quadrature (flat-facet, constant phase per face).
#   Upgrade to Gordon's exact triangle integral for coarse-mesh accuracy (Phase 4).
# ponytail: dot(ki,n)<0 illumination only — no ray-traced self-shadowing.
#   Add BVH occlusion query in Phase 5 for concave/re-entrant shapes.

References:
  Knott, Schaeffer, Tuley, *Radar Cross Section*, 2nd ed. (2004) §5.2
  Ufimtsev (2014), *Fundamentals of the Physical Theory of Diffraction*, §6
"""

using LinearAlgebra: norm
using KernelAbstractions

# ── Scalar helpers ────────────────────────────────────────────────────────────

@inline function _normalize3(v, ::Type{T}) where {T<:AbstractFloat}
    x, y, z = T(v[1]), T(v[2]), T(v[3])
    n = sqrt(x*x + y*y + z*z)
    n > eps(T) || throw(ArgumentError("zero-length direction vector"))
    return (x/n, y/n, z/n)
end

@inline _cross3(a::NTuple{3}, b::NTuple{3}) =
    (a[2]*b[3] - a[3]*b[2],
     a[3]*b[1] - a[1]*b[3],
     a[1]*b[2] - a[2]*b[1])

@inline _dot3(a::NTuple{3}, b::NTuple{3}) =
    a[1]*b[1] + a[2]*b[2] + a[3]*b[3]

# ── KernelAbstractions kernel ─────────────────────────────────────────────────

# Per-facet PO contributions.  Each thread handles one face i.
# Arguments:
#   Gre, Gim  — 3×Nf output (real/imag parts of complex G vectors)
#   normals   — 3×Nf unit outward normals (AbstractMatrix{T})
#   centroids — 3×Nf face centroids       (AbstractMatrix{T})
#   areas     — Nf face areas             (AbstractVector{T})
#   kice      — NTuple{3,T}: ki × ei (precomputed outside kernel)
#   ki        — NTuple{3,T}: incident direction (illumination test)
#   smki      — NTuple{3,T}: s − ki (phase direction)
#   k         — T: wavenumber
#
# Shadowed faces (dot(ki,n)≥0) write zero contributions.
# Shadowed faces contribute zero; sign: lit iff dot(ki,n) < 0 (normal toward source).
#
# ponytail: separate Gre/Gim arrays instead of Complex{T} for portability across GPU backends.
@kernel function _po_kernel!(
    Gre, Gim,
    normals, centroids, areas,
    kice, ki, smki, k,
)
    i  = @index(Global, Linear)
    nx = normals[1, i]; ny = normals[2, i]; nz = normals[3, i]

    if ki[1]*nx + ki[2]*ny + ki[3]*nz < 0
        # n_i × kice  (kice = ki × ei)
        cx  = ny*kice[3] - nz*kice[2]
        cy  = nz*kice[1] - nx*kice[3]
        cz  = nx*kice[2] - ny*kice[1]
        # phase = k · dot(s−ki, c_i)
        ph  = k * (smki[1]*centroids[1,i] + smki[2]*centroids[2,i] + smki[3]*centroids[3,i])
        sr  = cos(ph);  si = sin(ph)
        ai  = areas[i]
        Gre[1,i] = ai*cx*sr;  Gim[1,i] = ai*cx*si
        Gre[2,i] = ai*cy*sr;  Gim[2,i] = ai*cy*si
        Gre[3,i] = ai*cz*sr;  Gim[3,i] = ai*cz*si
    else
        # shadowed — zero contribution
        T = eltype(areas)
        z = zero(T)
        Gre[1,i] = z;  Gim[1,i] = z
        Gre[2,i] = z;  Gim[2,i] = z
        Gre[3,i] = z;  Gim[3,i] = z
    end
end

# ── Core dispatch ─────────────────────────────────────────────────────────────

"""
    _po_gperp_core(...) → (p1, p2, p3)::NTuple{3,Complex{T}}

Internal: run the PO KA kernel and return the complex transverse amplitude G_⊥
as a 3-tuple of Complex{T}.  Used by both `_po_rcs_core` and `_ptd_rcs_core`.
"""
function _po_gperp_core(
    normals   :: AbstractMatrix{T},
    centroids :: AbstractMatrix{T},
    areas     :: AbstractVector{T},
    ki  :: NTuple{3,T},
    s   :: NTuple{3,T},
    ei  :: NTuple{3,T},
    k   :: T,
) where {T<:AbstractFloat}
    Nf      = length(areas)
    backend = get_backend(normals)
    Gre     = similar(normals, T, (3, Nf))
    Gim     = similar(normals, T, (3, Nf))

    kice = _cross3(ki, ei)
    smki = (s[1]-ki[1], s[2]-ki[2], s[3]-ki[3])

    kern = _po_kernel!(backend, 256)
    kern(Gre, Gim, normals, centroids, areas, kice, ki, smki, k; ndrange = Nf)
    KernelAbstractions.synchronize(backend)

    # Reduce rows: sum over faces (axis 2)
    # ponytail: simple sequential reduction; use parallel tree-reduce on GPU in Phase 4.
    Gr1 = sum(view(Gre, 1, :));  Gr2 = sum(view(Gre, 2, :));  Gr3 = sum(view(Gre, 3, :))
    Gi1 = sum(view(Gim, 1, :));  Gi2 = sum(view(Gim, 2, :));  Gi3 = sum(view(Gim, 3, :))

    # Gperp = G − dot(s, G) · s   (s is real)
    dre = s[1]*Gr1 + s[2]*Gr2 + s[3]*Gr3
    dim = s[1]*Gi1 + s[2]*Gi2 + s[3]*Gi3

    p1 = complex(Gr1 - dre*s[1], Gi1 - dim*s[1])
    p2 = complex(Gr2 - dre*s[2], Gi2 - dim*s[2])
    p3 = complex(Gr3 - dre*s[3], Gi3 - dim*s[3])

    return (p1, p2, p3)
end

"""
    _po_rcs_core(normals, centroids, areas, ki, s, ei, k) → σ

Internal: run the KA kernel, reduce, compute σ.
All arguments are already normalised, converted to type T.
Arrays are passed raw so they can be CPU or GPU (Phase-4 ready).
"""
function _po_rcs_core(
    normals   :: AbstractMatrix{T},
    centroids :: AbstractMatrix{T},
    areas     :: AbstractVector{T},
    ki  :: NTuple{3,T},
    s   :: NTuple{3,T},
    ei  :: NTuple{3,T},
    k   :: T,
) where {T<:AbstractFloat}
    p1, p2, p3 = _po_gperp_core(normals, centroids, areas, ki, s, ei, k)
    return (k^2 / T(π)) * (abs2(p1) + abs2(p2) + abs2(p3))
end

# ── Public API ────────────────────────────────────────────────────────────────

"""
    po_rcs(mesh, ki, s, ei; k) → σ [m²]

Bistatic Physical Optics RCS of PEC mesh.

- `ki` — incident direction unit vector (direction of propagation)
- `s`  — scattering (observation) direction unit vector
- `ei` — incident E-field polarisation unit vector (must be ⊥ ki)
- `k`  — wavenumber 2π/λ [rad/m]

Returns monostatic RCS σ [m²].  Generic over float type of `mesh`.

# ponytail: TypeContracts not applied — standalone function, same decision as
#   Phase-1 scsn_general / flat_disk_rcs (contracts.jl comment).
"""
function po_rcs(
    mesh :: TriMesh{T},
    ki, s, ei;
    k    :: Real,
) where {T<:AbstractFloat}
    kiT = _normalize3(ki, T)
    sT  = _normalize3(s,  T)
    eiT = _normalize3(ei, T)
    abs(_dot3(kiT, eiT)) < T(0.1) ||
        throw(ArgumentError("ei is not ⊥ ki: |dot(ki,ei)| = $(abs(_dot3(kiT, eiT)))"))
    return _po_rcs_core(mesh.normals, mesh.centroids, mesh.areas, kiT, sT, eiT, T(k))
end

"""
    po_rcs_monostatic(mesh, ki, ei; k) → σ [m²]

Monostatic (backscatter) Physical Optics RCS of PEC mesh.

`s = −ki` is set automatically.  See `po_rcs` for argument descriptions.
"""
function po_rcs_monostatic(
    mesh :: TriMesh{T},
    ki, ei;
    k    :: Real,
) where {T<:AbstractFloat}
    kiT = _normalize3(ki, T)
    eiT = _normalize3(ei, T)
    abs(_dot3(kiT, eiT)) < T(0.1) ||
        throw(ArgumentError("ei is not ⊥ ki: |dot(ki,ei)| = $(abs(_dot3(kiT, eiT)))"))
    sT = (.-kiT[1], .-kiT[2], .-kiT[3])
    return _po_rcs_core(mesh.normals, mesh.centroids, mesh.areas, kiT, sT, eiT, T(k))
end

"""
    po_rcs_sweep(mesh, dirs, ei; k) → Vector{T}

Monostatic PO RCS over a collection of incidence directions.

- `dirs` — 3×N matrix; column j is the j-th incidence direction (normalised internally).
- `ei`   — polarisation unit vector (applied to all angles; warn if not ⊥ any dir).
- `k`    — wavenumber [rad/m].

Returns a `Vector{T}` of length N.  Mesh arrays are extracted once for efficiency.

# ponytail: serial loop over directions; parallelise across dirs in Phase 4 if
#   the sweep dominates (N >> 1000).
"""
function po_rcs_sweep(
    mesh :: TriMesh{T},
    dirs :: AbstractMatrix{<:Real},
    ei;
    k    :: Real,
) where {T<:AbstractFloat}
    size(dirs, 1) == 3 || throw(ArgumentError("dirs must be 3×N (got $(size(dirs,1))×…)"))
    eiT = _normalize3(ei, T)
    kT  = T(k)
    N   = size(dirs, 2)
    σs  = Vector{T}(undef, N)
    nrm, cen, ar = mesh.normals, mesh.centroids, mesh.areas
    for j in 1:N
        kiT  = _normalize3(view(dirs, :, j), T)
        sT   = (.-kiT[1], .-kiT[2], .-kiT[3])
        σs[j] = _po_rcs_core(nrm, cen, ar, kiT, sT, eiT, kT)
    end
    return σs
end

"""
    to_dbsm(σ) → σ_dBsm

Convert RCS from m² to dBsm (decibels relative to 1 m²).
"""
to_dbsm(σ::Real) = 10 * log10(σ)
