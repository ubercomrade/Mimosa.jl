using Test
using Mimosa

const SITEGA_FIXTURES = joinpath(@__DIR__, "..", "fixtures")

@testset "SiteGA constructor" begin
    # Basic construction: (25, motif_length) matrix
    rep = Matrix{Float32}(undef, 25, 10)
    fill!(rep, 0.0f0)
    m = SiteGA("test", rep, 10)
    @test m.name == "test"
    @test m.motif_length == 10
    @test size(m.representation) == (25, 10)
    @test length(m) == 10
    @test eltype(m) == Float32

    # Invalid: wrong row count
    @test_throws MimosaError SiteGA("bad", Matrix{Float32}(undef, 24, 3), 3)
    # Invalid: zero motif length
    @test_throws MimosaError SiteGA("bad", Matrix{Float32}(undef, 25, 0), 0)
    # Invalid: non-finite values
    bad_rep = Matrix{Float32}(undef, 25, 3)
    fill!(bad_rep, 0.0f0)
    bad_rep[1, 1] = NaN32
    @test_throws MimosaError SiteGA("bad", bad_rep, 3)
end

@testset "SiteGA show" begin
    rep = Matrix{Float32}(undef, 25, 12)
    fill!(rep, 0.0f0)
    m = SiteGA("gata2", rep, 12)
    s = sprint(show, m)
    @test contains(s, "SiteGA")
    @test contains(s, "gata2")
end

@testset "SiteGA equality" begin
    rep_a = Matrix{Float32}(undef, 25, 3)
    fill!(rep_a, 1.0f0)
    rep_b = Matrix{Float32}(undef, 25, 3)
    fill!(rep_b, 1.0f0)

    a = SiteGA("x", rep_a, 3)
    b = SiteGA("x", rep_b, 3)
    c = SiteGA("y", rep_a, 3)

    @test a == b
    @test a != c  # different name
    @test isapprox(a, b)
end

@testset "SiteGA scorebounds" begin
    rep = Matrix{Float32}(undef, 25, 3)
    fill!(rep, 0.0f0)
    # Set some non-trivial values
    rep[1, :] = [1.0f0, 2.0f0, 3.0f0]    # nuc1=0, nuc2=0 (aa)
    rep[7, :] = [-1.0f0, -2.0f0, -3.0f0]  # nuc1=1, nuc2=1 (cc)
    rep[13, :] = [0.5f0, 1.5f0, 2.5f0]   # nuc1=2, nuc2=2 (gg)
    m = SiteGA("test", rep, 3)

    mn, mx = scorebounds(m)
    # col1: min(-1, 0, 0.5, 0...) = -1.0, max = 1.0
    # col2: min(-2, 0, 1.5, 0...) = -2.0, max = 2.0
    # col3: min(-3, 0, 2.5, 0...) = -3.0, max = 3.0
    @test mn ≈ -6.0f0
    @test mx ≈ 6.0f0
end

@testset "SiteGA parsing" begin
    # Test reading sitega_gata2.mat
    path = joinpath(SITEGA_FIXTURES, "sitega_gata2.mat")
    @test isfile(path)

    m = read_sitega(path)
    @test m.name == "GATA2"
    @test m.motif_length == 12
    @test size(m.representation) == (25, 12)

    # Test reading sitega.mat
    path2 = joinpath(SITEGA_FIXTURES, "sitega.mat")
    @test isfile(path2)
    m2 = read_sitega(path2)
    @test m2.motif_length == 12

    # Test reading sitega_stat6.mat
    path3 = joinpath(SITEGA_FIXTURES, "sitega_stat6.mat")
    @test isfile(path3)
    m3 = read_sitega(path3)
    @test size(m3.representation) == (25, m3.motif_length)

    # File not found
    @test_throws MimosaError read_sitega("nonexistent.mat")
end

@testset "SiteGA scanning contract" begin
    model = read_sitega(joinpath(SITEGA_FIXTURES, "sitega_gata2.mat"))
    sequence = UInt8[0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3]
    _test_scan_contract(model, sequence)
    _test_batch_scan_contract(model, (sequence, reverse(sequence)))
end

@testset "SiteGA write round-trip" begin
    path = joinpath(SITEGA_FIXTURES, "sitega_gata2.mat")
    m = read_sitega(path)

    # Write to temp file and re-read
    tmp = tempname() * ".mat"
    write_sitega(tmp, m)
    m2 = read_sitega(tmp)

    @test m2.name == m.name
    @test m2.motif_length == m.motif_length
    @test m2.representation ≈ m.representation

    rm(tmp)
end
