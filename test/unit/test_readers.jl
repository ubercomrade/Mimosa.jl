using Test
using Mimosa

@testset "Low-level PFM readers are internal" begin
    @test !(:read_meme in names(Mimosa))
    @test !(:read_pfm in names(Mimosa))
end

const REPO_ROOT = dirname(dirname(@__DIR__))
const EXAMPLES = joinpath(REPO_ROOT, "examples")

@testset "read_meme parses pif4.meme" begin
    pfm = Mimosa.read_meme(joinpath(EXAMPLES, "pif4.meme"); index=0)
    @test pfm.name == "pwm_model"
    @test size(pfm.frequencies) == (4, 12)
    @test eltype(pfm.frequencies) === Float32
end

@testset "read_pfm parses pif4.pfm" begin
    pfm = Mimosa.read_pfm(joinpath(EXAMPLES, "pif4.pfm"))
    @test size(pfm.frequencies) == (4, 12)
    @test eltype(pfm.frequencies) === Float32
end

@testset "PWM reader adapters" begin
    meme_pwm = Mimosa._read_meme_pwm(joinpath(EXAMPLES, "pif4.meme"))
    pfm_pwm = Mimosa._read_pfm_pwm(joinpath(EXAMPLES, "pif4.pfm"))
    @test meme_pwm isa PWM
    @test pfm_pwm isa PWM
    @test meme_pwm.representation ≈ pfm_pwm.representation
end

@testset "readmodel detects format from extension" begin
    pwm_meme = readmodel(joinpath(EXAMPLES, "pif4.meme"))
    @test pwm_meme isa PWM
    @test size(pwm_meme.representation) == (5, 12)
    pwm_pfm = readmodel(joinpath(EXAMPLES, "pif4.pfm"))
    @test pwm_pfm isa PWM
    @test size(pwm_pfm.representation) == (5, 12)
    # Both should produce identical PWM weights since pif4.pfm and pif4.meme
    # represent the same motif.
    @test pwm_meme.representation ≈ pwm_pfm.representation
end

@testset "MEME malformed inputs" begin
    tmp = tempname()
    open(tmp, "w") do io
        println(io, "MEME version 4")
        println(io)
        println(io, "MOTIF")
    end
    @test_throws Mimosa.ModelFormatError Mimosa.read_meme(tmp)

    tmp2 = tempname()
    open(tmp2, "w") do io
        println(io, "MEME version 4")
        println(io)
        println(io, "MOTIF x")
        println(io, "letter-probability matrix: alength= 4 w= 3 nsites= 10")
        println(io, "0.25 0.25 0.25 0.25")
        println(io, "0.25 0.25 0.25")  # wrong column count
    end
    @test_throws Mimosa.ModelFormatError Mimosa.read_meme(tmp2)

    @test_throws Mimosa.ModelFormatError Mimosa.read_meme(
        joinpath(EXAMPLES, "pif4.meme"); index=99
    )
    @test_throws Mimosa.ModelFormatError Mimosa.read_meme("/nonexistent/path.meme")
end

@testset "PFM malformed inputs" begin
    tmp = tempname()
    open(tmp, "w") do io
        println(io, "0.1 0.2 0.3")
    end
    @test_throws Mimosa.ModelFormatError Mimosa.read_pfm(tmp)

    @test_throws Mimosa.ModelFormatError Mimosa.read_pfm("/nonexistent/path.pfm")
end

@testset "readmodel unsupported format" begin
    tmp = tempname() * ".txt"
    open(tmp, "w") do io
        println(io, "garbage")
    end
    @test_throws Mimosa.ModelFormatError readmodel(tmp)
end
