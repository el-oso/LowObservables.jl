"""
Physical Theory of Diffraction (PTD) edge correction — Phase 5b-ii.

Integrates Ufimtsev's elementary edge wave (EEW) directivity along all mesh edges
and adds it coherently to the Physical Optics far-field amplitude.

Convention: e^{+jωt} (engineering), same as physical_optics.jl.
The Ufimtsev EEW spatial phase (originally e^{-iωt}) is converted by i→−i:
  e^{ik(k̂ⁱ−ŝ)·r}  →  e^{ik(ŝ−k̂ⁱ)·r} = exp(im·k·dot(s−ki, r))  [Julia form]
This matches PO's exp(im·k·dot(s−ki, c_i)) exactly. ✓

Total far-field amplitude (coefficient of e^{ikR}/R):
  f_total = f_PO + Σ_edges f_edge
where:
  f_PO   = (im·k / 2π) · G_⊥^PO               [m, complex 3-vector]
  f_edge = (1/2π)·(E₀ₜ·F⁽¹⁾ + Z₀H₀ₜ·G⁽¹⁾)·L·sinc·e^{ik(s−ki)·c_e}  [m]
and σ = 4π|f_total|².

References: Ufimtsev (2014) §7.8, Eqs. 7.136–7.143.

Edges handled: 1-face open rims (α=2π half-plane) and 2-face convex dihedrals
(α=π+acos(n̂_A·n̂_B); cube→3π/2, fold→2π).  A dihedral is treated as a REAL edge
only if it is sharp enough (|α−π| ≥ `sharp_threshold`); near-coplanar tessellation
edges of a smooth body are excluded — see `_ptd_rcs_core` for why.

# ponytail: serial CPU edge loop; GPU edge kernel deferred.
# ponytail: concave edges (α<π from normal-dot formula) skipped; flag if needed.
# ponytail: non-manifold edges (>2 faces) skipped.
"""

"""
    ptd_rcs_monostatic(mesh, ki, ei; k) → σ [m²]

Monostatic PO+PTD RCS of a PEC triangulated surface.

Physical Optics field plus Ufimtsev elementary edge-wave contributions integrated
along all mesh edges.  Arguments identical to `po_rcs_monostatic`.

- `ki` — incident direction unit vector (propagation direction)
- `ei` — incident E-field polarisation unit vector (must be ⊥ ki)
- `k`  — wavenumber 2π/λ [rad/m]

Generic over float type of `mesh` (Float32 and Float64 supported).
"""
function ptd_rcs_monostatic(
    mesh :: TriMesh{T},
    ki, ei;
    k    :: Real,
    sharp_threshold :: Real = deg2rad(20),
) where {T<:AbstractFloat}
    kiT = _normalize3(ki, T)
    eiT = _normalize3(ei, T)
    abs(_dot3(kiT, eiT)) < T(0.1) ||
        throw(ArgumentError("ei is not ⊥ ki: |dot(ki,ei)| = $(abs(_dot3(kiT, eiT)))"))
    sT  = (.-kiT[1], .-kiT[2], .-kiT[3])   # monostatic: s = −ki
    _ptd_rcs_core(mesh, kiT, sT, eiT, T(k), T(sharp_threshold))
end

# Internal: PO + edge-wave coherent sum → σ.
#
# `sharp_threshold` [rad] — dihedral feature-edge criterion.  A 2-face edge is a
# REAL geometric edge (gets PTD fringe) only if |α−π| ≥ sharp_threshold; below it
# the faces are treated as a smooth surface (tessellation edge → no fringe).
# This is standard PTD-on-mesh practice: a coherent EEW sum over the facet edges
# of a smooth body does NOT cancel at finite mesh resolution (it stays ~O(1)),
# so near-coplanar edges MUST be excluded by a sharpness threshold rather than
# relied upon to self-cancel.  1-face (open rim) edges are always real (α=2π).
# ponytail: fixed 20° default; tune per target (lower for fine meshes / shallow
# real edges, higher to suppress coarse faceting).
function _ptd_rcs_core(
    mesh :: TriMesh{T},
    ki   :: NTuple{3,T},
    s    :: NTuple{3,T},
    ei   :: NTuple{3,T},
    k    :: T,
    sharp_threshold :: T,
) where {T<:AbstractFloat}

    # ── PO amplitude: f_PO = (im·k/2π) · G_⊥^PO ─────────────────────────────
    G1, G2, G3 = _po_gperp_core(mesh.normals, mesh.centroids, mesh.areas, ki, s, ei, k)
    fac = complex(zero(T), k / (2 * T(π)))   # = im·k/(2π)
    f1 = fac * G1
    f2 = fac * G2
    f3 = fac * G3

    # ── Edge corrections ──────────────────────────────────────────────────────
    verts = mesh.vertices

    for ((u, v), face_ids) in mesh.edge_faces
        nf = length(face_ids)
        nf == 0 && continue

        # ── Edge geometry ────────────────────────────────────────────────────
        ux = T(verts[1,u]); uy = T(verts[2,u]); uz = T(verts[3,u])
        vx = T(verts[1,v]); vy = T(verts[2,v]); vz = T(verts[3,v])
        ex = vx-ux; ey = vy-uy; ez = vz-uz
        L  = sqrt(ex*ex + ey*ey + ez*ez)
        L < eps(T) && continue
        cx_e = (ux+vx)/2; cy_e = (uy+vy)/2; cz_e = (uz+vz)/2   # midpoint
        t0x = ex/L; t0y = ey/L; t0z = ez/L  # raw edge tangent (sign fixed below)

        # ── Wedge angle α + lit face ─────────────────────────────────────────
        # Lit face = the one with dot(k̂ⁱ,n̂)<0.  Defines φ=0.  The other face (if
        # any) sits at φ=α.  Both lit → double-illumination (eew handles φ₀>π).
        local α::T
        local fi_lit::Int
        if nf == 1
            α      = 2 * T(π)                       # open rim → half-plane
            fi     = face_ids[1]
            (ki[1]*mesh.normals[1,fi] + ki[2]*mesh.normals[2,fi] + ki[3]*mesh.normals[3,fi] < zero(T)) || continue
            fi_lit = fi
        elseif nf == 2
            fi1 = face_ids[1]; fi2 = face_ids[2]
            n1x = T(mesh.normals[1,fi1]); n1y = T(mesh.normals[2,fi1]); n1z = T(mesh.normals[3,fi1])
            n2x = T(mesh.normals[1,fi2]); n2y = T(mesh.normals[2,fi2]); n2z = T(mesh.normals[3,fi2])
            d12 = clamp(n1x*n2x + n1y*n2y + n1z*n2z, T(-1), T(1))
            α   = T(π) + acos(d12)                  # convex dihedral (π→3π/2→2π)
            # ponytail: concave edge (α<π) skipped + flagged; needs face-ordering
            # to disambiguate inside/outside, deferred.
            α < T(π) - T(1e-6) && continue
            # Feature-edge test: near-coplanar faces = smooth surface, not a real
            # edge → no fringe (see _ptd_rcs_core header).  Also avoids the α=π pole.
            abs(α - T(π)) < sharp_threshold && continue
            lit1 = ki[1]*n1x + ki[2]*n1y + ki[3]*n1z < zero(T)
            lit2 = ki[1]*n2x + ki[2]*n2y + ki[3]*n2z < zero(T)
            (lit1 || lit2) || continue              # both shadowed → skip
            fi_lit = lit1 ? fi1 : fi2               # lit face defines φ=0
        else
            continue                                # non-manifold
        end

        nlx = T(mesh.normals[1,fi_lit]); nly = T(mesh.normals[2,fi_lit]); nlz = T(mesh.normals[3,fi_lit])

        # ── Edge frame (right-handed, ẑ=t̂, ŷ=n̂_lit on the vacuum side) ──────
        # x̂ = lit-face surface direction ⊥ t̂ pointing AWAY from edge (toward the
        # lit-face centroid) → this is φ=0.  ŷ = n̂_lit (vacuum side).  ẑ = x̂×ŷ
        # is the edge tangent with sign fixed for a right-handed (x̂,ŷ,ẑ) frame.
        dcx = T(mesh.centroids[1,fi_lit]) - cx_e
        dcy = T(mesh.centroids[2,fi_lit]) - cy_e
        dcz = T(mesh.centroids[3,fi_lit]) - cz_e
        dct = dcx*t0x + dcy*t0y + dcz*t0z
        xrx = dcx - dct*t0x; xry = dcy - dct*t0y; xrz = dcz - dct*t0z
        xn  = sqrt(xrx*xrx + xry*xry + xrz*xrz)
        xn < eps(T) && continue          # centroid on the edge line, degenerate
        xhx = xrx/xn; xhy = xry/xn; xhz = xrz/xn
        yhx = nlx; yhy = nly; yhz = nlz                       # ŷ = n̂_lit
        tx = xhy*yhz - xhz*yhy           # ẑ = t̂ = x̂ × ŷ
        ty = xhz*yhx - xhx*yhz
        tz = xhx*yhy - xhy*yhx

        # ── Cone angle γ₀ = acos(−dot(ki,t̂)) ───────────────────────────────
        g0     = acos(clamp(-(ki[1]*tx + ki[2]*ty + ki[3]*tz), T(-1), T(1)))
        sing0  = sin(g0)
        sing0 < T(1e-10) && continue     # ki ∥ edge, degenerate

        # ── Incidence azimuth φ₀ (source = −ki, projected ⊥ t̂) ─────────────
        mkix  = -ki[1]; mkiy = -ki[2]; mkiz = -ki[3]
        mkit  = mkix*tx + mkiy*ty + mkiz*tz       # dot(-ki, t̂)
        mkipx = mkix - mkit*tx; mkipy = mkiy - mkit*ty; mkipz = mkiz - mkit*tz
        mkipl = sqrt(mkipx*mkipx + mkipy*mkipy + mkipz*mkipz)
        mkipl < T(1e-10) && continue     # ki ∥ edge (caught above, but guard)
        dmkix = mkipx/mkipl; dmkiy = mkipy/mkipl; dmkiz = mkipz/mkipl
        phi0  = mod(atan(dmkix*yhx + dmkiy*yhy + dmkiz*yhz,
                         dmkix*xhx + dmkiy*xhy + dmkiz*xhz), T(2)*T(π))

        # ── Observation azimuth φ (ŝ projected ⊥ t̂) ────────────────────────
        st    = s[1]*tx + s[2]*ty + s[3]*tz         # dot(ŝ, t̂) = cos(ϑ)
        vartheta = acos(clamp(st, T(-1), T(1)))
        sinvth   = sin(vartheta)
        sinvth < T(1e-10) && continue    # ŝ ∥ edge, degenerate
        spx = s[1] - st*tx; spy = s[2] - st*ty; spz = s[3] - st*tz
        spl = sqrt(spx*spx + spy*spy + spz*spz)
        spl < T(1e-10) && continue
        dsx = spx/spl; dsy = spy/spl; dsz = spz/spl
        phi = mod(atan(dsx*yhx + dsy*yhy + dsz*yhz,
                       dsx*xhx + dsy*xhy + dsz*xhz), T(2)*T(π))

        # Validate azimuth ranges; skip if outside valid domain
        phi0 < T(1e-8) && (phi0 = T(1e-8))   # clamp off φ₀=0 singularity
        phi  < T(1e-8) && (phi  = T(1e-8))
        phi  > α - T(1e-8) && (phi = α - T(1e-8))
        phi0 > α           && continue         # source in shadow region — skip

        # ── EEW directivity ──────────────────────────────────────────────────
        fg = eew_directivity(phi, phi0, g0, vartheta, α)
        Ft = fg.F_theta   # Complex{T} or T (on-cone real)
        Gt = fg.G_theta
        Gp = fg.G_phi
        # Guard: fringe formula blows up at α=π (coplanar faces: cot singularity
        # at reflection boundary).  Physically fringe = 0 for flat junctions, so
        # skip any edge where EEW returns non-finite values.
        # ponytail: also catches Float32 cosh overflow in off-cone _sigma12 branch.
        (isfinite(real(Ft)) && isfinite(real(Gt)) && isfinite(real(Gp))) || continue

        # ── Tangential incident-field components (E₀ = 1) ───────────────────
        E0t = ei[1]*tx + ei[2]*ty + ei[3]*tz        # ê_i · t̂
        # k̂ⁱ × ê_i tangential: Z₀·H₀ₜ = (k̂ⁱ × ê_i) · t̂
        kicex = ki[2]*ei[3] - ki[3]*ei[2]
        kicey = ki[3]*ei[1] - ki[1]*ei[3]
        kicez = ki[1]*ei[2] - ki[2]*ei[1]
        H0t   = kicex*tx + kicey*ty + kicez*tz      # Z₀·H₀ₜ

        # ── Sinc factor: ∫_{−L/2}^{L/2} e^{ik(ŝ−k̂ⁱ)·t̂·ζ} dζ = L·sinc(u) ──
        smki_t = (s[1]-ki[1])*tx + (s[2]-ki[2])*ty + (s[3]-ki[3])*tz
        u      = k * smki_t * L / 2
        sinc_v = abs(u) < T(1e-8) ? one(T) : sin(u) / u

        # ── Propagation phase: exp(im·k·dot(ŝ−k̂ⁱ, c_e)) ────────────────────
        ph    = k * ((s[1]-ki[1])*cx_e + (s[2]-ki[2])*cy_e + (s[3]-ki[3])*cz_e)
        phase = complex(cos(ph), sin(ph))

        # ── Per-edge complex amplitude (scalar × Cartesian basis) ─────────────
        prefac    = L * sinc_v / (2 * T(π))           # real scalar [m/rad]
        # E0t and H0t are real T; Ft,Gt,Gp may be Complex{T} (off-cone) or T (on-cone)
        amp_theta = (E0t * Ft + H0t * Gt) * prefac * phase   # Complex{T}
        amp_phi   = (H0t * Gp)            * prefac * phase   # Complex{T}

        # ϑ̂_edge = cos(ϑ)·ρ̂ − sin(ϑ)·t̂,   φ̂_edge = t̂ × ρ̂
        thhx = st*dsx - sinvth*tx
        thhy = st*dsy - sinvth*ty
        thhz = st*dsz - sinvth*tz
        phx  = ty*dsz - tz*dsy     # t̂ × ρ̂
        phy  = tz*dsx - tx*dsz
        phz  = tx*dsy - ty*dsx

        # Accumulate into total far-field amplitude
        f1 += amp_theta * thhx + amp_phi * phx
        f2 += amp_theta * thhy + amp_phi * phy
        f3 += amp_theta * thhz + amp_phi * phz
    end

    # σ = 4π |f_total|²
    return 4 * T(π) * (abs2(f1) + abs2(f2) + abs2(f3))
end
