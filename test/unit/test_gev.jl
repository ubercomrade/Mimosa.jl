using Test
using Mimosa
using Random

@testset "GEV fit" begin
    # Basic fit on synthetic Gumbel data
    rng = MersenneTwister(42)
    gumbel_scores = [randn(rng) for _ in 1:200]  # Approx Gumbel-like
    result = fit_gev(gumbel_scores)

    @test result isa GEVFit
    @test result.scale > 0
    @test isfinite(result.shape)
    @test isfinite(result.location)
    @test isfinite(result.scale)

    # Shape should be near 0 for approximately Gumbel data
    @test abs(result.shape) < 0.5

    # Survival function properties
    @test 0.0 ≤ survival(result, minimum(gumbel_scores) - 100) ≤ 1.0
    @test 0.0 ≤ survival(result, maximum(gumbel_scores) + 100) ≤ 1.0
    @test survival(result, maximum(gumbel_scores) + 100) ≈ 0.0 atol = 1e-10
    @test survival(result, minimum(gumbel_scores) - 100) ≈ 1.0 atol = 1e-10

    # CDF + SF = 1
    for x in [-2.0, 0.0, 2.0]
        @test cdf(result, x) + survival(result, x) ≈ 1.0 atol = 1e-10
    end

    # scipy_params sign flip: c = -k
    c, loc, scale = scipy_params(result)
    @test c ≈ -result.shape
    @test loc ≈ result.location
    @test scale ≈ result.scale
end

@testset "GEV input validation" begin
    @test_throws ArgumentError fit_gev([1.0, 2.0, 3.0]; max_iter=0)
    @test_throws ArgumentError fit_gev([1.0, 2.0, 3.0]; tol=0.0)
    @test_throws ArgumentError GEVFit(0.0, 0.0, 0.0, true, 1, 0.0)
    @test_throws ArgumentError GEVFit(0.0, 0.0, 1.0, true, -1, 0.0)
    @test_throws ArgumentError GEVFit(NaN, 0.0, 1.0, true, 1, 0.0)
end

@testset "GEV fit edge cases" begin
    # Too few scores
    result = fit_gev([1.0, 2.0])
    @test result isa GEVFitFailure
    @test occursin("at least 3", result.message)

    # Constant sample
    result = fit_gev(fill(3.14, 10))
    @test result isa GEVFitFailure
    @test occursin("identical", result.message)

    # NaN in sample
    result = fit_gev([1.0, NaN, 3.0, 4.0, 5.0])
    @test result isa GEVFitFailure
    @test occursin("non-finite", result.message)

    # Inf in sample
    result = fit_gev([1.0, Inf, 3.0, 4.0, 5.0])
    @test result isa GEVFitFailure
    @test occursin("non-finite", result.message)

    # Empty sample
    result = fit_gev(Float64[])
    @test result isa GEVFitFailure
    @test occursin("at least 3", result.message)
end

@testset "GEV survival function edge cases" begin
    # Gumbel case (shape = 0)
    gev = GEVFit(0.0, 0.0, 1.0, true, 1, 0.0)
    @test survival(gev, 0.0) ≈ 1.0 - exp(-exp(0.0)) atol = 1e-10
    @test cdf(gev, 0.0) ≈ exp(-exp(0.0)) atol = 1e-10

    # Fréchet case (shape > 0): lower bounded
    gev = GEVFit(0.5, 0.0, 1.0, true, 1, 0.0)
    # Below lower bound: SF = 1
    @test survival(gev, -3.0) ≈ 1.0
    # At a normal point
    @test 0.0 < survival(gev, 1.0) < 1.0

    # Weibull case (shape < 0): upper bounded
    gev = GEVFit(-0.5, 0.0, 1.0, true, 1, 0.0)
    # Above upper bound: SF = 0
    @test survival(gev, 3.0) ≈ 0.0
    # At a normal point
    @test 0.0 < survival(gev, -1.0) < 1.0
end
