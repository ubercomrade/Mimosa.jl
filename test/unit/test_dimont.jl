using Test
using Mimosa

const DIMONT_FIXTURES = joinpath(@__DIR__, "..", "fixtures")

@testset "Dimont constructor" begin
    # span=0: (5, motif_length) matrix
    rep0 = Matrix{Float32}(undef, 5, 10)
    fill!(rep0, 0.0f0)
    m0 = Dimont("test0", rep0, 0, 10)
    @test m0.name == "test0"
    @test m0.span == 0
    @test m0.motif_length == 10
    @test size(m0.representation) == (5, 10)
    @test length(m0) == 10
    @test eltype(m0) == Float32

    # span=1: (25, motif_length) matrix
    rep1 = Matrix{Float32}(undef, 25, 8)
    fill!(rep1, 0.0f0)
    m1 = Dimont("test1", rep1, 1, 8)
    @test m1.span == 1
    @test size(m1.representation) == (25, 8)

    # span=3: (625, motif_length) matrix
    rep3 = Matrix{Float32}(undef, 625, 5)
    fill!(rep3, 0.0f0)
    m3 = Dimont("test3", rep3, 3, 5)
    @test m3.span == 3
    @test size(m3.representation) == (625, 5)

    # Invalid: wrong row count for span
    @test_throws MimosaError Dimont("bad", Matrix{Float32}(undef, 4, 3), 0, 3)
    @test_throws MimosaError Dimont("bad", Matrix{Float32}(undef, 24, 3), 1, 3)

    # Invalid: zero motif length
    @test_throws MimosaError Dimont("bad", Matrix{Float32}(undef, 5, 0), 0, 0)

    # Invalid: negative span
    @test_throws MimosaError Dimont("bad", Matrix{Float32}(undef, 5, 3), -1, 3)

    # Invalid: non-finite values
    bad_rep = Matrix{Float32}(undef, 5, 3)
    fill!(bad_rep, 0.0f0)
    bad_rep[1, 1] = NaN32
    @test_throws MimosaError Dimont("bad", bad_rep, 0, 3)
end

@testset "Dimont show" begin
    rep = Matrix{Float32}(undef, 5, 13)
    fill!(rep, 0.0f0)
    m = Dimont("example", rep, 0, 13)
    s = sprint(show, m)
    @test contains(s, "Dimont")
    @test contains(s, "example")
    @test contains(s, "span=0")
end

@testset "Dimont equality" begin
    rep_a = Matrix{Float32}(undef, 5, 3)
    fill!(rep_a, 1.0f0)
    rep_b = Matrix{Float32}(undef, 5, 3)
    fill!(rep_b, 1.0f0)

    a = Dimont("x", rep_a, 0, 3)
    b = Dimont("x", rep_b, 0, 3)
    c = Dimont("y", rep_a, 0, 3)
    rep_d = Matrix{Float32}(undef, 25, 3)
    fill!(rep_d, 1.0f0)
    d = Dimont("x", rep_d, 1, 3)

    @test a == b
    @test a != c  # different name
    @test a != d  # different span
    @test isapprox(a, b)
end

@testset "Dimont scorebounds" begin
    rep = Matrix{Float32}(undef, 5, 3)
    fill!(rep, 0.0f0)
    rep[1, :] = [1.0f0, 2.0f0, 3.0f0]    # A
    rep[2, :] = [-1.0f0, -2.0f0, -3.0f0]  # C
    rep[3, :] = [0.5f0, 1.5f0, 2.5f0]   # G
    rep[4, :] = [0.0f0, 0.0f0, 0.0f0]   # T
    rep[5, :] = [-1.0f0, -2.0f0, -3.0f0]  # N = min
    m = Dimont("test", rep, 0, 3)

    mn, mx = scorebounds(m)
    @test mn ≈ -6.0f0  # sum of col mins: -1 + -2 + -3
    @test mx ≈ 6.0f0   # sum of col maxes: 1 + 2 + 3
end

@testset "Dimont parsing" begin
    # Test reading exampleD-model-1.xml (span=0, length=13)
    path = joinpath(DIMONT_FIXTURES, "exampleD-model-1.xml")
    @test isfile(path)
    m = read_dimont(path)
    @test m.name == "exampleD-model-1"
    @test m.span == 0
    @test m.motif_length == 13
    @test size(m.representation) == (5, 13)

    # Test reading stat_dimont-model-1.xml (span=3, length=5)
    path2 = joinpath(DIMONT_FIXTURES, "stat_dimont-model-1.xml")
    @test isfile(path2)
    m2 = read_dimont(path2)
    @test m2.span == 3
    @test m2.motif_length == 5
    @test size(m2.representation) == (625, 5)

    # Test reading PEAKS036274 (span=3, length=10)
    path3 = joinpath(DIMONT_FIXTURES, "PEAKS036274_FOXA1_P35582_MACS2-model-1.xml")
    @test isfile(path3)
    m3 = read_dimont(path3)
    @test m3.span == 3
    @test m3.motif_length == 10
    @test size(m3.representation) == (625, 10)

    # File not found
    @test_throws MimosaError read_dimont("nonexistent.xml")
end

@testset "Dimont scanning single sequence" begin
    # Use stat_dimont (span=3, length=5, window=8, kmer=4)
    path = joinpath(DIMONT_FIXTURES, "stat_dimont-model-1.xml")
    m = read_dimont(path)

    seq = UInt8[0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3]
    n_pos = npositions(m, length(seq))
    # window = motif_length + span = 5 + 3 = 8
    @test n_pos == 20 - 8 + 1  # 13

    fwd = scan(m, seq; strands=ForwardOnly())
    @test length(fwd) == n_pos
    @test all(isfinite, fwd)

    rev = scan(m, seq; strands=ReverseOnly())
    @test length(rev) == n_pos
    @test all(isfinite, rev)

    best = scan(m, seq; strands=BestStrand())
    @test length(best) == n_pos
    @test all(best .>= min.(fwd, rev))

    both = scan(m, seq; strands=BothStrands())
    @test both.forward ≈ fwd
    @test both.reverse ≈ rev

    # In-place scan
    dest = Vector{Float32}(undef, n_pos)
    scan!(dest, m, seq; strands=ForwardOnly())
    @test dest ≈ fwd

    # Short sequence
    short_seq = UInt8[0, 1, 2, 3]
    fwd_short = scan(m, short_seq; strands=ForwardOnly())
    @test isempty(fwd_short)
end

@testset "Dimont scanning batch" begin
    path = joinpath(DIMONT_FIXTURES, "stat_dimont-model-1.xml")
    m = read_dimont(path)

    seqs = [
        UInt8[0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3],
        UInt8[3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0],
        UInt8[0, 0, 0, 0],
    ]
    batch = EncodedSequenceBatch(seqs)

    scores = scan(m, batch; strands=ForwardOnly())
    @test nrows(scores) == 3
    @test rowlength(scores, 1) == 20 - 8 + 1  # window=8
    @test rowlength(scores, 3) == 0  # short sequence

    both = scan(m, batch; strands=BothStrands())
    @test nrows(both.forward) == 3
    @test nrows(both.reverse) == 3

    best = scan(m, batch; strands=BestStrand())
    @test nrows(best) == 3

    lens = Mimosa.scan_result_lengths(m, batch)
    @test length(lens) == 3
    @test lens[1] == 20 - 8 + 1
    @test lens[3] == 0
end

@testset "Dimont scan determinism" begin
    path = joinpath(DIMONT_FIXTURES, "stat_dimont-model-1.xml")
    m = read_dimont(path)
    seq = UInt8[0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3]

    s1 = scan(m, seq; strands=ForwardOnly())
    s2 = scan(m, seq; strands=ForwardOnly())
    @test s1 == s2

    # Non-mutation
    seq_copy = copy(seq)
    scan(m, seq; strands=ForwardOnly())
    @test seq == seq_copy
end

@testset "Dimont span=0 equivalence to order-0 BaMM scan" begin
    # A span=0 Dimont should produce the same scanning geometry as an order=0 BaMM:
    # kmer=1, context=0, window=motif_length, n_terms=motif_length
    path = joinpath(DIMONT_FIXTURES, "exampleD-model-1.xml")
    m = read_dimont(path)

    @test Mimosa.kmer(m) == 1
    @test Mimosa.context_length(m) == 0
    @test Mimosa.window_size(m) == m.motif_length
    @test Mimosa.scan_width(m) == m.motif_length

    # Verify scanning produces finite results
    seq = UInt8[0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3]
    fwd = scan(m, seq; strands=ForwardOnly())
    @test length(fwd) == 20 - 13 + 1  # 8
    @test all(isfinite, fwd)
end
