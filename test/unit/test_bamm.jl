using Test
using Mimosa

const BAMM_FIXTURES = joinpath(@__DIR__, "..", "fixtures")

@testset "BaMM constructor" begin
    # Basic construction
    rep = Matrix{Float32}(undef, 5, 4)
    fill!(rep, 0.0f0)
    m = BaMM("test", rep, 0, 4)
    @test m.name == "test"
    @test m.order == 0
    @test m.motif_length == 4
    @test size(m.representation) == (5, 4)
    @test length(m) == 4
    @test eltype(m) == Float32

    # Order 1: 5^2 = 25 rows
    rep1 = Matrix{Float32}(undef, 25, 3)
    fill!(rep1, 0.0f0)
    m1 = BaMM("test1", rep1, 1, 3)
    @test size(m1.representation) == (25, 3)
    @test Mimosa.kmer(m1) == 2
    @test Mimosa.context_length(m1) == 1
    @test Mimosa.window_size(m1) == 4

    # Invalid: wrong row count
    @test_throws MimosaError BaMM("bad", Matrix{Float32}(undef, 4, 3), 0, 3)
    # Invalid: negative order
    @test_throws MimosaError BaMM("bad", Matrix{Float32}(undef, 5, 3), -1, 3)
    # Invalid: zero motif length
    @test_throws MimosaError BaMM("bad", Matrix{Float32}(undef, 5, 0), 0, 0)
    # Invalid: non-finite values
    bad_rep = Matrix{Float32}(undef, 5, 3)
    fill!(bad_rep, 0.0f0)
    bad_rep[1, 1] = NaN32
    @test_throws MimosaError BaMM("bad", bad_rep, 0, 3)
end

@testset "BaMM show" begin
    rep = Matrix{Float32}(undef, 5, 4)
    fill!(rep, 0.0f0)
    m = BaMM("myog", rep, 0, 4)
    s = sprint(show, m)
    @test contains(s, "BaMM")
    @test contains(s, "myog")
    @test contains(s, "order=0")
end

@testset "BaMM equality" begin
    rep_a = Matrix{Float32}(undef, 5, 3)
    fill!(rep_a, 1.0f0)
    rep_b = Matrix{Float32}(undef, 5, 3)
    fill!(rep_b, 1.0f0)

    a = BaMM("x", rep_a, 0, 3)
    b = BaMM("x", rep_b, 0, 3)
    c = BaMM("y", rep_a, 0, 3)

    # Different order requires different row count
    rep_d = Matrix{Float32}(undef, 25, 3)
    fill!(rep_d, 1.0f0)
    d = BaMM("x", rep_d, 1, 3)

    @test a == b
    @test a != c  # different name
    @test a != d  # different order

    @test isapprox(a, b)
end

@testset "BaMM scorebounds" begin
    rep = Matrix{Float32}(undef, 5, 3)
    rep[1, :] = [1.0f0, 2.0f0, 3.0f0]
    rep[2, :] = [-1.0f0, -2.0f0, -3.0f0]
    rep[3, :] = [0.0f0, 0.0f0, 0.0f0]
    rep[4, :] = [0.5f0, 1.5f0, 2.5f0]
    rep[5, :] = [-2.0f0, -3.0f0, -4.0f0]  # N-state min
    m = BaMM("test", rep, 0, 3)

    mn, mx = scorebounds(m)
    @test mn ≈ -9.0f0  # sum of per-column min: -2 + -3 + -4
    # col1: max(1, -1, 0, 0.5, -2) = 1.0
    # col2: max(2, -2, 0, 1.5, -3) = 2.0
    # col3: max(3, -3, 0, 2.5, -4) = 3.0
    @test mx ≈ 6.0f0
end

@testset "BaMM parsing" begin
    # Test reading myog.ihbcp
    path = joinpath(BAMM_FIXTURES, "myog.ihbcp")
    @test isfile(path)

    # Default order (max order from file)
    m = read_bamm(path)
    @test m.name == "myog"
    @test m.motif_length == 14
    @test size(m.representation) == (5^(m.order + 1), 14)

    # Order 0
    m0 = read_bamm(path; order=0)
    @test m0.order == 0
    @test size(m0.representation) == (5, 14)

    # Order 1
    m1 = read_bamm(path; order=1)
    @test m1.order == 1
    @test size(m1.representation) == (25, 14)

    # Order exceeds max → clamped to max
    m_clamped = read_bamm(path; order=100)
    @test m_clamped.order == m.order

    # File not found
    @test_throws MimosaError read_bamm("nonexistent.ihbcp")
end

@testset "BaMM scanning single sequence" begin
    path = joinpath(BAMM_FIXTURES, "myog.ihbcp")
    m1 = read_bamm(path; order=1)

    # Create a simple test sequence (longer than window_size = 14 + 1 = 15)
    seq = UInt8[0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3]
    n_pos = npositions(m1, length(seq))
    @test n_pos == 20 - 15 + 1  # seq_len - window_size + 1

    # Forward scan
    fwd = scan(m1, seq; strands=ForwardOnly())
    @test length(fwd) == n_pos
    @test all(isfinite, fwd)

    # Reverse scan
    rev = scan(m1, seq; strands=ReverseOnly())
    @test length(rev) == n_pos
    @test all(isfinite, rev)

    # Best strand
    best = scan(m1, seq; strands=BestStrand())
    @test length(best) == n_pos
    @test all(best .>= min.(fwd, rev))

    # Both strands
    both = scan(m1, seq; strands=BothStrands())
    @test both.forward ≈ fwd
    @test both.reverse ≈ rev

    # In-place scan
    dest = Vector{Float32}(undef, n_pos)
    scan!(dest, m1, seq; strands=ForwardOnly())
    @test dest ≈ fwd

    # Short sequence (shorter than window_size)
    short_seq = UInt8[0, 1, 2, 3]
    fwd_short = scan(m1, short_seq; strands=ForwardOnly())
    @test isempty(fwd_short)
end

@testset "BaMM scanning batch" begin
    path = joinpath(BAMM_FIXTURES, "myog.ihbcp")
    m0 = read_bamm(path; order=0)
    m1 = read_bamm(path; order=1)

    # Build a batch
    seqs = [
        UInt8[0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3],
        UInt8[3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0],
        UInt8[0, 0, 0, 0],
    ]
    batch = EncodedSequenceBatch(seqs)

    # Forward scan order=0
    scores0 = scan(m0, batch; strands=ForwardOnly())
    @test nrows(scores0) == 3
    @test rowlength(scores0, 1) == 20 - 14 + 1  # order=0: window=14
    @test rowlength(scores0, 3) == 0  # short sequence

    # Forward scan order=1
    scores1 = scan(m1, batch; strands=ForwardOnly())
    @test nrows(scores1) == 3
    @test rowlength(scores1, 1) == 20 - 15 + 1  # order=1: window=15
    @test rowlength(scores1, 3) == 0  # short sequence

    # Both strands
    both = scan(m1, batch; strands=BothStrands())
    @test nrows(both.forward) == 3
    @test nrows(both.reverse) == 3

    # Best strand
    best = scan(m1, batch; strands=BestStrand())
    @test nrows(best) == 3

    # Scan result lengths
    lens = Mimosa.scan_result_lengths(m1, batch)
    @test length(lens) == 3
    @test lens[1] == 20 - 15 + 1
    @test lens[3] == 0
end

@testset "BaMM order=0 equivalence to PWM scan" begin
    # When order=0, BaMM scanning should behave like PWM scanning
    # (each position scored independently, kmer=1, context=0)
    path = joinpath(BAMM_FIXTURES, "myog.ihbcp")
    m0 = read_bamm(path; order=0)

    seq = UInt8[0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3]

    # BaMM order=0 scan
    bamm_scores = scan(m0, seq; strands=ForwardOnly())

    # Manual PWM-style scan using the 5-row representation
    W = m0.motif_length
    n_pos = max(length(seq) - W + 1, 0)
    pwm_scores = Vector{Float32}(undef, n_pos)
    for pos in 1:n_pos
        total = 0.0f0
        for p in 1:W
            base = Int(seq[pos + p - 1]) + 1
            total += m0.representation[base, p]
        end
        pwm_scores[pos] = total
    end

    @test bamm_scores ≈ pwm_scores
end

@testset "BaMM scan determinism" begin
    path = joinpath(BAMM_FIXTURES, "myog.ihbcp")
    m1 = read_bamm(path; order=1)
    seq = UInt8[0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3]

    s1 = scan(m1, seq; strands=ForwardOnly())
    s2 = scan(m1, seq; strands=ForwardOnly())
    @test s1 == s2

    # Non-mutation: scanning should not modify the model or sequence
    seq_copy = copy(seq)
    scan(m1, seq; strands=ForwardOnly())
    @test seq == seq_copy
end
