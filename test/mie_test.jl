# Mie series validation for PEC sphere RCS.
#
# Reference values from Bohren & Huffman (1983) and the well-known Mie curve
# for a PEC sphere. See also: Wiscombe (1980) for series truncation.

@testitem "Rayleigh limit: σ ≪ πa² for ka ≪ 1" begin
    using LowObservables
    a = 0.5   # 50 cm sphere
    for ka in [0.01, 0.05, 0.1]
        λ = 2π * a / ka
        σ = mie_rcs_sphere(a, λ)
        σ_opt = optical_limit_rcs(a)
        # In the Rayleigh region σ/(πa²) ≈ 9(ka)^4 → very small
        @test σ / σ_opt < 1e-2
    end
end

@testitem "Rayleigh limit: matches 9(ka)^4 to 0.5%" begin
    using LowObservables
    a = 1.0
    # The leading Rayleigh term for a PEC sphere: σ/(πa²) = 9(ka)^4
    for ka in [0.01, 0.02, 0.05]
        λ = 2π * a / ka
        σ = mie_rcs_sphere(a, λ)
        σ_opt = optical_limit_rcs(a)
        predicted = 9 * ka^4
        @test isapprox(σ / σ_opt, predicted; rtol = 5e-3)
    end
end

@testitem "First Mie resonance peak near ka≈1" begin
    using LowObservables
    # The σ/(πa²) curve for a PEC sphere has its first maximum near ka≈1.03,
    # with a normalized value of ≈3.65 (Bohren & Huffman §4.4).
    a = 1.0
    # ka=1.0 is in the peak region; test that we're in the known [3.0, 5.0] band
    λ = 2π * a / 1.0
    σ = mie_rcs_sphere(a, λ)
    ratio = σ / optical_limit_rcs(a)
    @test 3.0 ≤ ratio ≤ 5.0

    # Tighter: the first peak value is ~3.65, test within 10%
    λ_peak = 2π * a / 1.03
    σ_peak = mie_rcs_sphere(a, λ_peak)
    ratio_peak = σ_peak / optical_limit_rcs(a)
    @test isapprox(ratio_peak, 3.655; rtol = 0.10)
end

@testitem "Optical limit: window average σ/(πa²) ≈ 1 within 5%" begin
    using LowObservables
    # Average over ka ∈ [30, 31] to smooth out Mie ripple.
    # The mean should converge to 1.0 (geometric optics).
    a = 0.3
    kas = range(30.0, 31.0, length=200)
    σ_values = [mie_rcs_sphere(a, 2π * a / ka) for ka in kas]
    σ_opt = optical_limit_rcs(a)
    mean_ratio = sum(σ_values) / (length(σ_values) * σ_opt)
    @test isapprox(mean_ratio, 1.0; atol = 0.05)
end

@testitem "Float type genericity: Float32 in → Float32 out" begin
    using LowObservables
    a32 = Float32(0.1)
    λ32 = Float32(0.03)
    σ = mie_rcs_sphere(a32, λ32)
    @test σ isa Float32
    # Result should be consistent with Float64 to within Float32 precision
    σ64 = mie_rcs_sphere(Float64(a32), Float64(λ32))
    @test isapprox(Float64(σ), σ64; rtol = 1e-5)
end

@testitem "rcs_vs_size returns a vector of the same length" begin
    using LowObservables
    a = 0.1
    λs = collect(range(0.01, 0.1, length=50))
    σs = rcs_vs_size(a, λs)
    @test length(σs) == 50
    @test all(>(0), σs)
end

@testitem "optical_limit_rcs returns πa²" begin
    using LowObservables
    for a in [0.1, 0.5, 1.0, 2.5]
        @test isapprox(optical_limit_rcs(a), π * a^2)
    end
end

@testitem "MieSphere type satisfies AbstractRCSModel contract" begin
    using LowObservables, TypeContracts
    @test implements(MieSphere, AbstractRCSModel)
end

@testitem "MieSphere rcs and mie_rcs_sphere are consistent" begin
    using LowObservables
    m = MieSphere(0.2)
    λ = 0.05
    @test isapprox(rcs(m, λ), mie_rcs_sphere(0.2, λ))
end

@testitem "wiscombe_nmax gives sensible truncation counts" begin
    using LowObservables
    # For small ka, minimum is 5
    @test wiscombe_nmax(0.01) == 5
    # For larger ka, scales up
    @test wiscombe_nmax(100.0) > 100
    # Non-negative always
    @test wiscombe_nmax(1e-6) ≥ 5
end
