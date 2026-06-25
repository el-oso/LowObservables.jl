"""
TypeContracts.jl interface definitions for LowObservables.

Users implement `AbstractRCSModel` to add custom RCS models.
"""

"""
    AbstractRCSModel

Interface contract for radar cross-section models.

Mandatory methods:
- `rcs(model, λ) :: Real` — monostatic RCS [m²] at wavelength λ [m]

Optional methods:
- `rcs_vs_wavelength(model, λs) :: AbstractVector` — RCS over a vector of wavelengths
"""
abstract type AbstractRCSModel end

@contract AbstractRCSModel begin
    rcs(::Self, λ::Real)::Real
    :optional
    rcs_vs_wavelength(::Self, λs::AbstractVector{<:Real})::AbstractVector{<:Real}
end

"""
    rcs(model::AbstractRCSModel, λ::Real) -> Real

Monostatic (backscatter) radar cross-section [m²] at wavelength λ [m].
"""
function rcs end
