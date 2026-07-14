using Test
using Mimosa

const EXAMPLES = joinpath(dirname(dirname(@__DIR__)), "..", "examples")

@testset "PWM construction from frequency matrices" begin
    frequencies = Float32[0.25 0.5; 0.25 0.25; 0.25 0.1; 0.25 0.15]
    pwm = pwm_from_pfm(frequencies; background=0.25, name="test")
    @test length(pwm) == 2
    @test size(pwm) == (5, 2)
    @test eltype(pwm) === Float32
    @test pwm.name == "test"
    @test pwm.background === (0.25f0, 0.25f0, 0.25f0, 0.25f0)
end

@testset "PWM weights validation" begin
    # 6 rows — wrong nucleotide axis
    bad = Float32[1 2; 3 4; 5 6; 7 8; 9 10; 11 12]
    @test_throws Mimosa.ModelDimensionError PWM("x", bad, (0.25f0, 0.25f0, 0.25f0, 0.25f0))

    # 4 rows — missing N row
    four = Float32[1 2; 3 4; 5 6; 7 8]
    @test_throws Mimosa.ModelDimensionError PWM("x", four, (0.25f0, 0.25f0, 0.25f0, 0.25f0))
end

@testset "pfm_to_pwm and pcm_to_pfm" begin
    pfm = Float32[0.0 1.0; 0.0 0.0; 0.0 0.0; 1.0 0.0]
    pwm = pfm_to_pwm(pfm)
    @test size(pwm) == (4, 2)
    # log((1 + 1e-4) / 0.25) at [1,1]=0 and [4,2]=0; log((0 + 1e-4)/0.25) elsewhere
    @test pwm[1, 1] ≈ log((0.0f0 + Float32(1e-4)) / 0.25f0)
    @test pwm[4, 1] ≈ log((1.0f0 + Float32(1e-4)) / 0.25f0)

    pcm = Float32[10 0; 0 0; 0 0; 0 10]
    pfm2 = pcm_to_pfm(pcm)
    @test sum(pfm2; dims=1) ≈ ones(Float32, 1, 2)
end

@testset "reverse_complement involutive" begin
    frequencies = Float32[0.1 0.2 0.3; 0.4 0.5 0.6; 0.7 0.8 0.9; 1.0 0.0 0.5]
    rc = reverse_complement(frequencies)
    @test reverse_complement(rc) ≈ frequencies
    @test size(rc) == (4, 3)
    # A row maps to T row reversed
    @test rc[4, :] ≈ reverse(frequencies[1, :])
end

@testset "scorebounds" begin
    pwm = PWM("t", Float32[1 2; 3 4; 5 6; 7 8; -100 -200], (0.25f0, 0.25f0, 0.25f0, 0.25f0))
    mn, mx = scorebounds(pwm)
    @test mn == sum([-100, -200])
    @test mx == sum([7, 8])
end
