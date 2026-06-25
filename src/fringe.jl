"""
2-D fringe directivity functions f⁽¹⁾ and g⁽¹⁾ for first-order PTD.

These are the *finite* directivity coefficients of the edge-diffracted (fringe) wave:
  f⁽¹⁾ — soft/E-polarisation (Dirichlet BC)
  g⁽¹⁾ — hard/H-polarisation (Neumann BC)

They are constructed so that the PO singularities at the GO boundaries (shadow boundary
ψ₁=π and reflection boundary ψ₂=π) cancel, leaving finite values everywhere.

Source: Ufimtsev (2014), *Fundamentals of the Physical Theory of Diffraction*, 2nd ed.,
MATLAB listing `fun_fg.m`, pp. 447–448.

# ponytail: TypeContracts @contract skipped — standalone function, same rationale as
# bodies_of_revolution.jl.  Add if an AbstractFringeModel type hierarchy emerges.

# ponytail: fringe *integral* Int_calcFringe.m (surface integration of Eq 4.19, needed
# for far-field amplitude) is deferred — not required for directivity-pattern validation.
"""

"""
    fringe_directivity(φ, φ₀, α) -> (; f1, g1)

Fringe directivity functions f⁽¹⁾ (`f1`, soft/E-pol) and g⁽¹⁾ (`g1`, hard/H-pol)
for a perfectly-conducting wedge of exterior angle `α`.

Faithful port of `fun_fg.m`, Ufimtsev (2014), pp. 447–448.

# Arguments
- `φ`  : observation angle [rad], in `[0, α]`
- `φ₀` : incidence angle [rad], in `[ε, π−ε]` (measured from the illuminated face)
- `α`  : exterior wedge angle [rad]; `α > π` for a convex edge (e.g. `3π/2` for 270°)

# Returns NamedTuple
- `f1` : soft (Dirichlet) fringe directivity — finite at all GO boundaries
- `g1` : hard (Neumann) fringe directivity — finite at all GO boundaries

Both `f1` and `g1` are **finite** at the geometric-optics boundaries `φ = π − φ₀`
(reflection) and `φ = π + φ₀` (shadow): the PO divergences cancel in `f = f_total − f_PO`.

Generic over float type: `Float32` inputs → `Float32` outputs.

# Branch structure (matches MATLAB exactly)
- Single-side illumination: `ε ≤ φ₀ ≤ α − π − ε`  (requires `α > π`)
  - Grazing sub-branches for `ψ₁ ≈ π` and `ψ₂ ≈ π`
- Double-side illumination: `α − π + ε ≤ φ₀ ≤ π − ε`
  - Grazing sub-branches for `ψ₂ ≈ π` and `2α − ψ₂ ≈ π`
"""
function fringe_directivity(φ::Real, φ₀::Real, α::Real)
    T = promote_type(typeof(float(φ)), typeof(float(φ₀)), typeof(float(α)))
    _fringe_directivity(T(φ), T(φ₀), T(α))
end

function _fringe_directivity(φ::T, φ₀::T, α::T) where {T<:AbstractFloat}
    ε  = T(1e-5)               # eps from fun_fg.m: 0 < ε << 1 rad
    ψ₁ = φ  - φ₀
    ψ₂ = φ  + φ₀
    n  = α  / T(π)
    a1 = (T(π) - ψ₂) / (2n)
    a2 = (T(π) - ψ₁) / (2n)
    a3 = (T(π) + ψ₂) / (2n)
    a4 = (T(π) + ψ₁) / (2n)

    if φ₀ >= ε && φ₀ <= α - T(π) - ε           # ── single-side illumination ───

        if ψ₁ <= T(π) + ε && ψ₁ >= T(π) - ε     # shadow boundary ψ₁ ≈ π (angle2 → 0)
            # cot(a2) and cot((π−ψ₁)/2) cancel; remaining terms finite
            f1 =  (1/(2n)) * (cot(a3) - cot(a4)) - T(0.5) * cot((T(π) - ψ₂) / 2)
            g1 = -(1/(2n)) * (cot(a3) + cot(a4)) + T(0.5) * cot((T(π) - ψ₂) / 2)

        elseif ψ₂ <= T(π) + ε && ψ₂ >= T(π) - ε  # reflection boundary ψ₂ ≈ π (angle1 → 0)
            # cot(a1) and cot((π−ψ₂)/2) cancel; remaining terms finite
            f1 =  (1/(2n)) * (-cot(a2) + cot(a3) - cot(a4)) + T(0.5) * cot((T(π) - ψ₁) / 2)
            g1 = -(1/(2n)) * (cot(a2) + cot(a3) + cot(a4)) + T(0.5) * cot((T(π) - ψ₁) / 2)

        else                                        # general single-side
            f  =  (1/(2n)) * (cot(a1) - cot(a2) + cot(a3) - cot(a4))
            g  = -(1/(2n)) * (cot(a1) + cot(a2) + cot(a3) + cot(a4))
            f0 =  T(0.5) * (cot((T(π) - ψ₂) / 2) - cot((T(π) - ψ₁) / 2))
            g0 = -T(0.5) * (cot((T(π) - ψ₂) / 2) + cot((T(π) - ψ₁) / 2))
            f1 = f - f0
            g1 = g - g0
        end

    elseif φ₀ >= α - T(π) + ε && φ₀ <= T(π) - ε  # ── double-side illumination ──

        if ψ₂ <= T(π) + ε && ψ₂ >= T(π) - ε       # reflection boundary ψ₂ ≈ π (angle1 → 0)
            f1 =  (1/(2n)) * (-cot(a2) + cot(a3) - cot(a4)) -
                  T(0.5) * cot((T(π) - 2α + ψ₂) / 2)
            g1 = -(1/(2n)) * (cot(a2) + cot(a3) + cot(a4)) +
                  T(0.5) * cot((T(π) - 2α + ψ₂) / 2)

        elseif 2α - ψ₂ <= T(π) + ε && 2α - ψ₂ >= T(π) - ε  # lower-face boundary (angle3 → π)
            f1 =  (1/(2n)) * (cot(a1) - cot(a2) - cot(a4)) - T(0.5) * cot((T(π) - ψ₂) / 2)
            g1 = -(1/(2n)) * (cot(a1) + cot(a2) + cot(a4)) + T(0.5) * cot((T(π) - ψ₂) / 2)

        else                                        # general double-side
            f  =  (1/(2n)) * (cot(a1) - cot(a2) + cot(a3) - cot(a4))
            g  = -(1/(2n)) * (cot(a1) + cot(a2) + cot(a3) + cot(a4))
            f0 =  T(0.5) * (cot((T(π) - ψ₂) / 2) + cot((T(π) - 2α + ψ₂) / 2))
            g0 = -T(0.5) * (cot((T(π) - ψ₂) / 2) + cot((T(π) - 2α + ψ₂) / 2))
            f1 = f - f0
            g1 = g - g0
        end

    elseif φ₀ == α && φ == T(π) / 2             # ── grazing on upper face ────
        f1 =  (1/(2n)) * (-cot(a2) + cot(a3) - cot(a4))
        g1 = -(1/(2n)) * (cot(a2) + cot(a3) + cot(a4))

    else
        # Out-of-range: φ₀ near α−π boundary between single/double illumination,
        # or other edge cases. Return zeros (not a valid diffraction geometry).
        f1 = zero(T)
        g1 = zero(T)
    end

    return (; f1, g1)
end
