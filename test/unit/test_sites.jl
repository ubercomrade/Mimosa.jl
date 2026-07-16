using Test
using Mimosa

@testset "SiteHit" begin
    h = SiteHit(1, 5, Int8(0), 3.5f0)
    @test h.seq_index == 1
    @test h.start == 5
    @test h.strand == 0
    @test h.score == 3.5f0
end

@testset "SiteCollection" begin
    coll = SiteCollection([1, 2], [5, 10], Int8[0, 1], [3.0f0, 5.0f0])
    @test length(coll) == 2
    @test !isempty(coll)

    empty = empty_site_collection()
    @test isempty(empty)
    @test length(empty) == 0

    # Mismatched lengths should error
    @test_throws ArgumentError SiteCollection([1], [5, 10], Int8[0], [3.0f0])
end

@testset "BestPerSequence selection" begin
    # Build a simple PWM for testing
    pwm4 = Float32[0.1 0.2 0.3; 0.4 0.1 0.2; 0.3 0.5 0.1; 0.2 0.2 0.4]
    pwm = pwm_from_pfm(pwm4; name="test")

    # Build a batch with known sequences
    # Sequence 1: ACG (positions 1-3 match motif)
    # Sequence 2: TTT (reverse strand)
    seqs = [
        UInt8[0, 1, 2, 0, 1, 2, 0],  # ACGACGA
        UInt8[3, 3, 3, 3, 0, 0, 0],  # TTTTAAA (reverse complement of first 3 = TTT)
    ]
    batch = EncodedSequenceBatch(seqs)

    selector = BestPerSequence()
    coll = selectsites(pwm, batch, selector; strands=BothStrands())

    @test length(coll) == 2
    # Each sequence should have exactly one hit
    @test coll.seq_indices == [1, 2]
    # Scores should be finite
    @test all(isfinite, coll.scores)
end

@testset "ThresholdHits selection" begin
    pwm4 = Float32[0.1 0.2 0.3; 0.4 0.1 0.2; 0.3 0.5 0.1; 0.2 0.2 0.4]
    pwm = pwm_from_pfm(pwm4; name="test")

    seqs = [UInt8[0, 1, 2, 0, 1, 2, 0, 1, 2]]
    batch = EncodedSequenceBatch(seqs)

    # First scan to find the score range
    fwd = scan(pwm, batch; strands=ForwardOnly())
    max_score = maximum(fwd.data)

    # Threshold at max should give exactly 1 hit (or more if tied)
    coll = selectsites(pwm, batch, ThresholdHits(max_score); strands=ForwardOnly())
    @test length(coll) >= 1

    # Very low threshold should give all positions
    coll_all = selectsites(pwm, batch, ThresholdHits(Float32(-Inf)); strands=ForwardOnly())
    @test length(coll_all) == 7  # 9 - 3 + 1 = 7 positions

    # Both strands with BothStrands
    coll_both = selectsites(pwm, batch, ThresholdHits(Float32(-Inf)); strands=BothStrands())
    @test length(coll_both) == 14  # 7 forward + 7 reverse
end

@testset "sort_hits!" begin
    coll = SiteCollection(
        [3, 1, 2, 1], [10, 5, 15, 8], Int8[0, 1, 0, 0], [1.0f0, 3.0f0, 2.0f0, 3.0f0]
    )
    sort_hits!(coll)

    # Expected order: (seq=1, score=3), (seq=1, score=3 → but start=5 < 8, strand=1 > 0)
    # Actually: sort by (seq_index asc, -score asc, start asc, strand asc)
    # (1, -3, 5, 1) → first, (1, -3, 8, 0) → second, (2, -2, 15, 0) → third, (3, -1, 10, 0) → fourth
    @test coll.seq_indices == [1, 1, 2, 3]
    @test coll.scores == [3.0f0, 3.0f0, 2.0f0, 1.0f0]
    @test coll.starts == [5, 8, 15, 10]
    @test coll.strands == Int8[1, 0, 0, 0]
end

@testset "select_top_fraction" begin
    coll = SiteCollection(
        [1, 2, 3, 4, 5],
        [10, 20, 30, 40, 50],
        Int8[0, 0, 0, 0, 0],
        [1.0f0, 5.0f0, 3.0f0, 4.0f0, 2.0f0],
    )

    # Top 40% of 5 = max(1, floor(5*0.4)) = 2
    result = select_top_fraction(coll, 0.4)
    @test length(result) == 2
    # Top 2 by score: 5.0 and 4.0
    @test result.scores == [5.0f0, 4.0f0]

    # Top 100% = all
    result_all = select_top_fraction(coll, 1.0)
    @test length(result_all) == 5

    # Top 10% of 5 = max(1, floor(0.5)) = 1
    result_one = select_top_fraction(coll, 0.1)
    @test length(result_one) == 1
    @test result_one.scores[1] == 5.0f0
end

@testset "extract_site_matrix" begin
    # Build a batch
    seqs = [
        UInt8[0, 1, 2, 3, 0, 1, 2],  # ACGTACG
        UInt8[3, 2, 1, 0, 3, 2, 1],  # TGCATGC
    ]
    batch = EncodedSequenceBatch(seqs)

    # Forward strand hit at position 1 in sequence 1, width 3
    coll = SiteCollection([1], [1], Int8[0], [1.0f0])
    sites = extract_site_matrix(batch, coll, 3)
    @test sites[:, 1] == UInt8[0, 1, 2]  # ACG

    # Forward strand hit at position 2 in sequence 1, width 3
    coll2 = SiteCollection([1], [2], Int8[0], [1.0f0])
    sites2 = extract_site_matrix(batch, coll2, 3)
    @test sites2[:, 1] == UInt8[1, 2, 3]  # CGT

    # Reverse strand hit at position 1 in sequence 1, width 3
    # Forward site = ACG (0,1,2), reverse complement = CGT (1,2,3)
    coll3 = SiteCollection([1], [1], Int8[1], [1.0f0])
    sites3 = extract_site_matrix(batch, coll3, 3)
    # Reverse complement of [0,1,2] = [2,1,0] reversed and complemented
    # Original: A(0) C(1) G(2) → reverse: G(2) C(1) A(0) → complement: C(1) G(2) T(3)
    @test sites3[:, 1] == UInt8[1, 2, 3]  # CGT

    # Reverse strand hit at position 1 in sequence 2, width 3
    # Forward: T(3) G(2) C(1) → reverse: C(1) G(2) T(3) → complement: G(2) C(1) A(0)
    coll4 = SiteCollection([2], [1], Int8[1], [1.0f0])
    sites4 = extract_site_matrix(batch, coll4, 3)
    @test sites4[:, 1] == UInt8[2, 1, 0]  # GCA
end

@testset "build_pcm" begin
    # 4 sites of width 4: ACGT, ACGT, TGCA, TGCA
    # sites[p, h] = base at position p for hit h
    sites = UInt8[
        0 0 3 3
        1 1 2 2
        2 2 1 1
        3 3 0 0
    ]
    pcm = build_pcm(sites, 4)
    # Position 1 (row 1): [0, 0, 3, 3] → A=2, T=2
    @test pcm[1, 1] == 2.0f0  # A at position 1
    @test pcm[4, 1] == 2.0f0  # T at position 1
    # Position 2 (row 2): [1, 1, 2, 2] → C=2, G=2
    @test pcm[2, 2] == 2.0f0  # C at position 2
    @test pcm[3, 2] == 2.0f0  # G at position 2
    # Position 3 (row 3): [2, 2, 1, 1] → G=2, C=2
    @test pcm[3, 3] == 2.0f0
    @test pcm[2, 3] == 2.0f0
    # Position 4 (row 4): [3, 3, 0, 0] → T=2, A=2
    @test pcm[4, 4] == 2.0f0
    @test pcm[1, 4] == 2.0f0
    @test sum(pcm) == 16.0f0  # 4 sites × 4 positions
end

@testset "site_strings" begin
    sites = UInt8[
        0 3
        1 2
        2 1
        3 0
    ]
    strs = site_strings(sites)
    @test strs == ["ACGT", "TGCA"]

    # N handling: 1×2 matrix → 2 sites of width 1
    sites_n = UInt8[4 0]
    strs_n = site_strings(sites_n)
    @test strs_n == ["N", "A"]
end

@testset "reconstruct_pfm basic" begin
    # Build a simple PWM
    pwm4 = Float32[0.9 0.1 0.1; 0.05 0.8 0.1; 0.03 0.05 0.7; 0.02 0.05 0.1]
    pwm = pwm_from_pfm(pwm4; name="test")

    # Build sequences that should produce clear hits
    seqs = [UInt8[0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2]]  # ACG repeated
    batch = EncodedSequenceBatch(seqs)

    pfm = reconstruct_pfm(pwm, batch, BestPerSequence(); pseudocount=0.25f0)
    @test size(pfm) == (4, 3)
    # Each column should sum to 1.0
    col_sums = vec(sum(pfm; dims=1))
    @test all(isapprox.(col_sums, ones(Float32, 3); atol=Float32(1e-5)))
end

@testset "reconstruct_pfm composes keyword selectors" begin
    pwm = PWM("test", Float32[1 1; 0 0; 0 0; 0 0; 0 0], (0.25f0, 0.25f0, 0.25f0, 0.25f0))
    batch = EncodedSequenceBatch([UInt8[0, 0, 0]])

    by_keywords = reconstruct_pfm(
        pwm, batch; mode=:threshold, score_threshold=0.0f0, top_fraction=1.0
    )
    explicit = reconstruct_pfm(pwm, batch, TopFractionHits(1.0, ThresholdHits(0.0f0)))
    @test by_keywords == explicit
end

@testset "Serial/threaded site workflow equivalence" begin
    frequencies = Float32[0.1 0.2 0.3; 0.4 0.1 0.2; 0.3 0.5 0.1; 0.2 0.2 0.4]
    pwm = pwm_from_pfm(frequencies; name="threaded_sites")
    batch = EncodedSequenceBatch([
        encode_sequence("ACGTACGTACGT"),
        encode_sequence("TTTTGGGGCCCC"),
        encode_sequence("NNNNACGTNNNN"),
    ])

    for selector in (BestPerSequence(), ThresholdHits(Float32(-Inf)))
        serial = selectsites(
            pwm, batch, selector; strands=BothStrands(), execution=SerialExecution()
        )
        threaded = selectsites(
            pwm, batch, selector; strands=BothStrands(), execution=ThreadedExecution(4)
        )
        @test threaded == serial
    end

    serial_pfm = reconstruct_pfm(pwm, batch, BestPerSequence(); execution=SerialExecution())
    threaded_pfm = reconstruct_pfm(
        pwm, batch, BestPerSequence(); execution=ThreadedExecution(4)
    )
    @test threaded_pfm == serial_pfm
end

@testset "selectsites empty batch" begin
    pwm4 = Float32[0.25 0.25; 0.25 0.25; 0.25 0.25; 0.25 0.25]
    pwm = pwm_from_pfm(pwm4; name="test")

    empty_batch = empty_sequence_batch()
    coll = selectsites(pwm, empty_batch, BestPerSequence())
    @test isempty(coll)
end

@testset "selectsites short sequences" begin
    pwm4 = Float32[0.9 0.1 0.1 0.1; 0.05 0.8 0.1 0.1; 0.03 0.05 0.7 0.1; 0.02 0.05 0.1 0.7]
    pwm = pwm_from_pfm(pwm4; name="test")

    # Sequence shorter than motif width (4) → no hits
    seqs = [UInt8[0, 1, 2]]  # ACG, length 3 < 4
    batch = EncodedSequenceBatch(seqs)
    coll = selectsites(pwm, batch, BestPerSequence())
    @test isempty(coll)
end

@testset "selectsites ForwardOnly" begin
    pwm4 = Float32[0.9 0.1; 0.05 0.8; 0.03 0.05; 0.02 0.05]
    pwm = pwm_from_pfm(pwm4; name="test")

    seqs = [UInt8[0, 1, 0, 1, 0]]
    batch = EncodedSequenceBatch(seqs)

    coll_fwd = selectsites(pwm, batch, BestPerSequence(); strands=ForwardOnly())
    @test length(coll_fwd) == 1
    @test coll_fwd.strands[1] == 0  # forward
end

@testset "selectsites ReverseOnly" begin
    pwm4 = Float32[0.9 0.1; 0.05 0.8; 0.03 0.05; 0.02 0.05]
    pwm = pwm_from_pfm(pwm4; name="test")

    seqs = [UInt8[0, 1, 0, 1, 0]]
    batch = EncodedSequenceBatch(seqs)

    coll_rev = selectsites(pwm, batch, BestPerSequence(); strands=ReverseOnly())
    @test length(coll_rev) == 1
    @test coll_rev.strands[1] == 1  # reverse
end

@testset "Higher-order model sites (BaMM)" begin
    # Create a simple BaMM with order=1
    rep = Matrix{Float32}(undef, 25, 5)
    for i in 1:25
        for j in 1:5
            rep[i, j] = Float32(randn())
        end
    end
    model = BaMM("test_bamm", rep, 1, 5)
    @test site_start_offset(model) == 1
    @test length(model) == 5

    # Create a batch with sequences long enough for scanning
    # window_size = motif_length + order = 6, need seq_len >= 6
    batch = make_random_sequences(5, 20; seed=42)

    # Select sites (best per sequence)
    coll = selectsites(model, batch, BestPerSequence(); strands=BothStrands())
    @test length(coll) <= 5  # at most one per sequence
    @test length(coll) >= 0

    if length(coll) > 0
        # Verify site extraction with offset
        offset = site_start_offset(model)
        for i in 1:length(coll)
            seq = sequence(batch, coll.seq_indices[i])
            # The motif starts at scan_pos + offset
            motif_start = coll.starts[i] + offset
            # Motif window should fit within the sequence
            @test motif_start >= 1
            @test motif_start + length(model) - 1 <= length(seq)
        end

        # Reconstruct PFM
        pfm = reconstruct_pfm(model, batch, BestPerSequence(); pseudocount=0.1f0)
        @test size(pfm) == (4, 5)
        # Each column should sum to ~1 (with pseudocount)
        for col in 1:5
            @test sum(pfm[:, col]) ≈ 1.0f0 atol = 0.01
        end
    end
end

@testset "Higher-order model sites (SiteGA)" begin
    rep = Matrix{Float32}(undef, 25, 6)
    for i in 1:25
        for j in 1:6
            rep[i, j] = Float32(randn())
        end
    end
    model = SiteGA("test_sitega", rep, 6)
    @test site_start_offset(model) == 0
    @test length(model) == 6

    batch = make_random_sequences(5, 20; seed=42)

    coll = selectsites(model, batch, BestPerSequence(); strands=BothStrands())
    @test length(coll) <= 5

    if length(coll) > 0
        # SiteGA offset = 0, so motif starts at scan_pos directly
        for i in 1:length(coll)
            seq = sequence(batch, coll.seq_indices[i])
            @test coll.starts[i] >= 1
            @test coll.starts[i] + length(model) - 1 <= length(seq)
        end

        pfm = reconstruct_pfm(model, batch, BestPerSequence(); pseudocount=0.1f0)
        @test size(pfm) == (4, 6)
    end
end

@testset "Higher-order model sites (Dimont)" begin
    span = 1
    rep = Matrix{Float32}(undef, 5^(span + 1), 5)
    for i in 1:size(rep, 1)
        for j in 1:size(rep, 2)
            rep[i, j] = Float32(randn())
        end
    end
    model = Dimont("test_dimont", rep, span, 5)
    @test site_start_offset(model) == 1
    @test length(model) == 5

    batch = make_random_sequences(5, 20; seed=42)

    coll = selectsites(model, batch, BestPerSequence(); strands=BothStrands())
    @test length(coll) <= 5

    if length(coll) > 0
        offset = site_start_offset(model)
        for i in 1:length(coll)
            seq = sequence(batch, coll.seq_indices[i])
            motif_start = coll.starts[i] + offset
            @test motif_start >= 1
            @test motif_start + length(model) - 1 <= length(seq)
        end

        pfm = reconstruct_pfm(model, batch, BestPerSequence(); pseudocount=0.1f0)
        @test size(pfm) == (4, 5)
    end
end

@testset "Higher-order model sites (Slim)" begin
    span = 1
    rep = Matrix{Float32}(undef, 5^(span + 1), 5)
    for i in 1:size(rep, 1)
        for j in 1:size(rep, 2)
            rep[i, j] = Float32(randn())
        end
    end
    model = Slim("test_slim", rep, span, 5)
    @test site_start_offset(model) == 1
    @test length(model) == 5

    batch = make_random_sequences(5, 20; seed=42)

    coll = selectsites(model, batch, BestPerSequence(); strands=BothStrands())
    @test length(coll) <= 5

    if length(coll) > 0
        pfm = reconstruct_pfm(model, batch, BestPerSequence(); pseudocount=0.1f0)
        @test size(pfm) == (4, 5)
    end
end
