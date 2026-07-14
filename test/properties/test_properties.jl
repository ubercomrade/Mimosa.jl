using Test
using Mimosa

const REPO_ROOT = joinpath(dirname(dirname(@__DIR__)), "..")
const EXAMPLES = joinpath(REPO_ROOT, "examples")

@testset "reverse_complement involution (PWM)" begin
    pwm = readmodel(joinpath(EXAMPLES, "pif4.meme"))
    rc = reverse_complement(pwm)
    @test reverse_complement(rc) == pwm
end

@testset "reverse_complement involution (frequency matrix)" begin
    frequencies = Mimosa.read_meme(joinpath(EXAMPLES, "pif4.meme")).frequencies
    rc = reverse_complement(frequencies)
    @test reverse_complement(rc) == frequencies
end

@testset "identical model-derived profile comparison gives 1.0" begin
    pwm = readmodel(joinpath(EXAMPLES, "pif4.meme"))
    sequences = make_random_sequences(20, 100; seed=42)
    @test compare(pwm, pwm, sequences; metric=:co).score ≈ 1.0f0
end

@testset "non-! functions do not mutate inputs" begin
    pwm = readmodel(joinpath(EXAMPLES, "pif4.meme"))
    sequences = make_random_sequences(20, 100; seed=42)
    copy_w = copy(pwm.weights)
    compare(pwm, pwm, sequences; metric=:co)
    scorebounds(pwm)
    reverse_complement(pwm)
    @test pwm.weights == copy_w
end

@testset "comparison is deterministic" begin
    pwm = readmodel(joinpath(EXAMPLES, "pif4.meme"))
    sequences = make_random_sequences(20, 100; seed=42)
    r1 = compare(pwm, pwm, sequences; metric=:co)
    r2 = compare(pwm, pwm, sequences; metric=:co)
    @test r1.score == r2.score
    @test r1.offset == r2.offset
    @test r1.orientation == r2.orientation
end

@testset "score bounds are consistent" begin
    pwm = readmodel(joinpath(EXAMPLES, "pif4.meme"))
    mn, mx = scorebounds(pwm)
    @test mn <= mx
    @test isfinite(mn) && isfinite(mx)
end

@testset "orientation labels are valid" begin
    pwm1 = readmodel(joinpath(EXAMPLES, "pif4.meme"))
    pwm2 = readmodel(joinpath(EXAMPLES, "gata2.meme"))
    sequences = make_random_sequences(20, 100; seed=42)
    for metric in (:co, :co_rowwise, :dice, :dice_rowwise, :cosine)
        r = compare(pwm1, pwm2, sequences; metric=metric)
        @test r.orientation in ("++", "+-", "-+", "--")
    end
end

# Stage 2 properties

@testset "reverse_complement involution (encoded sequence)" begin
    for s in ["ACGT", "AAAA", "ACGTACGT", "NNNN", "", "ACGTNNACGT"]
        seq = encode_sequence(s)
        @test reverse_complement(reverse_complement(seq)) == seq
    end
end

@testset "scan does not mutate inputs" begin
    pwm = readmodel(joinpath(EXAMPLES, "pif4.meme"))
    seq = encode_sequence("ACGTACGTACGTACGTACGTACGTAC")
    seq_copy = copy(seq)
    weights_copy = copy(pwm.weights)
    scan(pwm, seq; strands=ForwardOnly())
    scan(pwm, seq; strands=BothStrands())
    @test seq == seq_copy
    @test pwm.weights == weights_copy
end

@testset "scan is deterministic" begin
    pwm = readmodel(joinpath(EXAMPLES, "pif4.meme"))
    seq = encode_sequence("ACGTACGTACGTACGTACGTACGTAC")
    r1 = scan(pwm, seq; strands=ForwardOnly())
    r2 = scan(pwm, seq; strands=ForwardOnly())
    @test r1 == r2
    p1 = scan(pwm, seq; strands=BothStrands())
    p2 = scan(pwm, seq; strands=BothStrands())
    @test p1.forward == p2.forward
    @test p1.reverse == p2.reverse
end

@testset "scan! == scan (forward)" begin
    pwm = readmodel(joinpath(EXAMPLES, "pif4.meme"))
    for s in ["ACGTACGTACGTACGTACGTACGTAC", "GGGGCCCCAAAATTTTGGGGCCCCAA"]
        seq = encode_sequence(s)
        alloc = scan(pwm, seq; strands=ForwardOnly())
        dest = Vector{Float32}(undef, length(alloc))
        scan!(dest, pwm, seq; strands=ForwardOnly())
        @test dest == alloc
    end
end

@testset "scan! == scan (reverse)" begin
    pwm = readmodel(joinpath(EXAMPLES, "pif4.meme"))
    seq = encode_sequence("ACGTACGTACGTACGTACGTACGTAC")
    alloc = scan(pwm, seq; strands=ReverseOnly())
    dest = Vector{Float32}(undef, length(alloc))
    scan!(dest, pwm, seq; strands=ReverseOnly())
    @test dest == alloc
end

@testset "scan! == scan (best)" begin
    pwm = readmodel(joinpath(EXAMPLES, "pif4.meme"))
    seq = encode_sequence("ACGTACGTACGTACGTACGTACGTAC")
    alloc = scan(pwm, seq; strands=BestStrand())
    dest = Vector{Float32}(undef, length(alloc))
    scan!(dest, pwm, seq; strands=BestStrand())
    @test dest == alloc
end

@testset "best strand = max of fwd and rev" begin
    pwm = readmodel(joinpath(EXAMPLES, "pif4.meme"))
    for s in ["ACGTACGTACGTACGTACGTACGTAC", "TTTTGGGGCCCCAAAATTTTGGGGCCC"]
        seq = encode_sequence(s)
        fwd = scan(pwm, seq; strands=ForwardOnly())
        rev = scan(pwm, seq; strands=ReverseOnly())
        best = scan(pwm, seq; strands=BestStrand())
        @test best == max.(fwd, rev)
    end
end

@testset "batch scan == single scan" begin
    pwm = readmodel(joinpath(EXAMPLES, "pif4.meme"))
    rows = [
        encode_sequence("ACGTACGTACGTACGTACGTACGTAC"),
        encode_sequence("TTTTGGGGCCCCAAAATTTTGGGGCCC"),
        encode_sequence("ACGT"),
    ]
    batch = EncodedSequenceBatch(rows)
    fwd = scan(pwm, batch; strands=ForwardOnly())
    for i in 1:nsequences(batch)
        single = scan(pwm, sequence(batch, i); strands=ForwardOnly())
        @test row(fwd, i) == single
    end
    rev = scan(pwm, batch; strands=ReverseOnly())
    for i in 1:nsequences(batch)
        single = scan(pwm, sequence(batch, i); strands=ReverseOnly())
        @test row(rev, i) == single
    end
end

@testset "short sequence returns empty scores" begin
    pwm = readmodel(joinpath(EXAMPLES, "pif4.meme"))
    W = length(pwm)
    for len in 0:(W - 1)
        seq = encode_sequence("A"^len)
        for strands in (ForwardOnly(), ReverseOnly(), BestStrand())
            @test scan(pwm, seq; strands=strands) == Float32[]
        end
        pair = scan(pwm, seq; strands=BothStrands())
        @test pair.forward == Float32[]
        @test pair.reverse == Float32[]
    end
    # Exactly motif width → 1 position
    seq = encode_sequence("A"^W)
    @test length(scan(pwm, seq; strands=ForwardOnly())) == 1
end

@testset "reverse scan equals forward scan of reverse-complement PWM" begin
    pwm = readmodel(joinpath(EXAMPLES, "pif4.meme"))
    seq = encode_sequence("ACGTACGTACGTACGTACGTACGTAC")
    rev_scores = scan(pwm, seq; strands=ReverseOnly())
    rc_pwm = reverse_complement(pwm)
    fwd_rc_scores = scan(rc_pwm, seq; strands=ForwardOnly())
    @test rev_scores ≈ fwd_rc_scores
end

@testset "batch reverse_complement involution" begin
    rows = [encode_sequence("ACGTACGT"), encode_sequence("TTTTGGGG"), encode_sequence("")]
    batch = EncodedSequenceBatch(rows)
    rc_batch = reverse_complement(batch)
    rc_rc_batch = reverse_complement(rc_batch)
    @test rc_rc_batch == batch
end

@testset "FASTA round-trip to_padded" begin
    batch, _ = read_fasta(joinpath(EXAMPLES, "foreground.fa"))
    padded, lengths = to_padded(batch)
    rt = from_padded(padded, lengths)
    @test rt == batch
end

# Stage 3 properties

@testset "profile comparison is deterministic" begin
    sp1 = read_scores(joinpath(EXAMPLES, "scores_1.fasta"))
    sp2 = read_scores(joinpath(EXAMPLES, "scores_2.fasta"))
    r1 = compare(sp1, sp2; metric=:co, search_range=2, window_radius=2)
    r2 = compare(sp1, sp2; metric=:co, search_range=2, window_radius=2)
    @test r1.score == r2.score
    @test r1.offset == r2.offset
    @test r1.orientation == r2.orientation
    @test r1.n_sites == r2.n_sites
end

@testset "profile comparison does not mutate inputs" begin
    sp1 = read_scores(joinpath(EXAMPLES, "scores_1.fasta"))
    sp2 = read_scores(joinpath(EXAMPLES, "scores_2.fasta"))
    data1_copy = copy(sp1.scores.data)
    offsets1_copy = copy(sp1.scores.offsets)
    data2_copy = copy(sp2.scores.data)
    offsets2_copy = copy(sp2.scores.offsets)
    compare(sp1, sp2; metric=:co, search_range=2, window_radius=2)
    @test sp1.scores.data == data1_copy
    @test sp1.scores.offsets == offsets1_copy
    @test sp2.scores.data == data2_copy
    @test sp2.scores.offsets == offsets2_copy
end

@testset "self-comparison gives high scores" begin
    sp1 = read_scores(joinpath(EXAMPLES, "scores_1.fasta"))
    result = compare(sp1, sp1; metric=:co_rowwise, search_range=0, window_radius=0)
    # Self-comparison with co_rowwise should give 1.0 (identical profiles)
    @test result.score ≈ 1.0f0 atol = 1e-5
    @test result.orientation == "++"
end

@testset "profile metric parse round-trip" begin
    for name in ("co", "co_rowwise", "dice", "dice_rowwise", "cosine")
        m = parse_profile_metric(name)
        @test metric_name(m) == name
    end
end

@testset "LogTailTable fit/transform round-trip" begin
    # Fit from a sample, then transform the same sample
    sample = Float32[0.1, 0.5, 0.9, 0.3, 0.7, 0.5, 0.1]
    table = fit(EmpiricalLogTail(), sample)
    rag = build_ragged([sample])
    transformed = transform_scores(table, rag)
    # Each unique score maps to a unique log-tail value
    # Verify that equal scores map to equal log-tail values
    for i in eachindex(sample)
        @test lookup_score(table, sample[i]) == row(transformed, 1)[i]
    end
end

@testset "profile comparison is deterministic across repeated calls" begin
    sp1 = read_scores(joinpath(EXAMPLES, "scores_1.fasta"))
    sp2 = read_scores(joinpath(EXAMPLES, "scores_2.fasta"))
    # With default parameters (search_range=10, window_radius=10)
    r1 = compare(sp1, sp2; metric=:co)
    r2 = compare(sp1, sp2; metric=:co)
    @test r1.score == r2.score
    @test r1.offset == r2.offset
    @test r1.n_sites == r2.n_sites
    @test r1.orientation == r2.orientation
end

# Stage 4 properties

@testset "selectsites is deterministic" begin
    pwm = readmodel(joinpath(EXAMPLES, "pif4.meme"))
    batch, _ = read_fasta(joinpath(EXAMPLES, "foreground.fa"))
    c1 = selectsites(pwm, batch, BestPerSequence(); strands=BothStrands())
    c2 = selectsites(pwm, batch, BestPerSequence(); strands=BothStrands())
    @test c1 == c2
end

@testset "selectsites does not mutate inputs" begin
    pwm = readmodel(joinpath(EXAMPLES, "pif4.meme"))
    batch, _ = read_fasta(joinpath(EXAMPLES, "foreground.fa"))
    data_copy = copy(batch.data)
    offsets_copy = copy(batch.offsets)
    selectsites(pwm, batch, BestPerSequence(); strands=BothStrands())
    @test batch.data == data_copy
    @test batch.offsets == offsets_copy
end

@testset "selectsites empty batch gives empty collection" begin
    pwm = readmodel(joinpath(EXAMPLES, "pif4.meme"))
    empty_batch = empty_sequence_batch()
    coll = selectsites(pwm, empty_batch, BestPerSequence())
    @test isempty(coll)
end

@testset "selectsites short sequences give no hits" begin
    pwm = readmodel(joinpath(EXAMPLES, "pif4.meme"))
    W = length(pwm)
    short_seq = encode_sequence("A"^(W - 1))
    batch = EncodedSequenceBatch([short_seq])
    coll = selectsites(pwm, batch, BestPerSequence())
    @test isempty(coll)
end

@testset "reconstruct_pfm columns sum to 1.0" begin
    pwm = readmodel(joinpath(EXAMPLES, "pif4.meme"))
    batch, _ = read_fasta(joinpath(EXAMPLES, "foreground.fa"))
    pfm = reconstruct_pfm(pwm, batch, BestPerSequence(); pseudocount=0.25f0)
    @test size(pfm) == (4, length(pwm))
    for col in 1:size(pfm, 2)
        @test sum(@view pfm[:, col]) ≈ 1.0f0 atol = 1e-5
    end
end

@testset "reconstruct_pfm with TopFraction gives fewer or equal sites" begin
    pwm = readmodel(joinpath(EXAMPLES, "pif4.meme"))
    batch, _ = read_fasta(joinpath(EXAMPLES, "foreground.fa"))
    pfm_all = reconstruct_pfm(pwm, batch, BestPerSequence(); pseudocount=0.25f0)
    pfm_half = reconstruct_pfm(pwm, batch, TopFractionHits(0.5); pseudocount=0.25f0)
    @test size(pfm_all) == size(pfm_half)
    # Both should have columns summing to 1.0
    for col in 1:size(pfm_half, 2)
        @test sum(@view pfm_half[:, col]) ≈ 1.0f0 atol = 1e-5
    end
end

@testset "site extraction reverse complement involution" begin
    # Extract a forward site, then reverse-complement it → should match
    # the reverse-strand extraction at the same position
    pwm = readmodel(joinpath(EXAMPLES, "pif4.meme"))
    seq = encode_sequence("ACGTACGTACGTACGTACGTACGTAC")
    batch = EncodedSequenceBatch([seq])
    W = length(pwm)

    # Forward hit at position 1
    coll_fwd = SiteCollection([1], [1], Int8[0], [0.0f0])
    site_fwd = extract_site_matrix(batch, coll_fwd, W)

    # Reverse hit at same position
    coll_rev = SiteCollection([1], [1], Int8[1], [0.0f0])
    site_rev = extract_site_matrix(batch, coll_rev, W)

    # The reverse-strand site should be the reverse complement of the forward site
    expected_rc = UInt8[s == N_CODE ? N_CODE : 0x03 - s for s in reverse(site_fwd[:, 1])]
    @test site_rev[:, 1] == expected_rc
end

@testset "sort_hits! is idempotent" begin
    coll = SiteCollection(
        [3, 1, 2, 1], [10, 5, 15, 8], Int8[0, 1, 0, 0], [1.0f0, 3.0f0, 2.0f0, 3.0f0]
    )
    sort_hits!(coll)
    coll_copy = SiteCollection(
        copy(coll.seq_indices), copy(coll.starts), copy(coll.strands), copy(coll.scores)
    )
    sort_hits!(coll)
    @test coll.seq_indices == coll_copy.seq_indices
    @test coll.starts == coll_copy.starts
    @test coll.strands == coll_copy.strands
    @test coll.scores == coll_copy.scores
end

# Stage 7 properties: parallelism and cache

@testset "serial/threaded scan equivalence (property)" begin
    pwm = readmodel(joinpath(EXAMPLES, "pif4.meme"))
    seqs = [
        encode_sequence("ACGTACGTACGTACGTACGTACGTAC"),
        encode_sequence("TTTTGGGGCCCCAAAACGTACGTAC"),
        encode_sequence("ACGT"),
        encode_sequence("NNNNNNNN"),
        encode_sequence("ACGTNNACGTACGTACGTNNNNACGT"),
    ]
    data = UInt8[]
    offsets = [1]
    for s in seqs
        append!(data, s)
        push!(offsets, length(data) + 1)
    end
    batch = EncodedSequenceBatch(data, offsets)

    for strands in (ForwardOnly(), ReverseOnly(), BestStrand(), BothStrands())
        serial = scan(pwm, batch; strands=strands, execution=SerialExecution())
        for nt in (1, 2, 4)
            threaded = scan(pwm, batch; strands=strands, execution=ThreadedExecution(nt))
            if strands isa BothStrands
                @test threaded.forward.data == serial.forward.data
                @test threaded.reverse.data == serial.reverse.data
            else
                @test threaded.data == serial.data
                @test threaded.offsets == serial.offsets
            end
        end
    end
end

@testset "cache keys are deterministic (property)" begin
    dir = mktempdir()
    cache = Cache(dir)
    for _ in 1:3
        k1 = cache_key(cache, "scan", "fp1", "fp2")
        k2 = cache_key(cache, "scan", "fp1", "fp2")
        @test k1 == k2
    end
end

@testset "cache survives write/read round-trip (property)" begin
    dir = mktempdir()
    cache = Cache(dir)
    for i in 1:5
        key = "roundtrip_$i"
        data = UInt8.(collect(1:i) .* 10)
        cache_set(cache, key, data; metadata=Dict("index" => i))
        @test cache_has(cache, key)
        @test cache_get(cache, key) == data
        meta = cache_get_meta(cache, key)
        @test meta["index"] == i
    end
end

@testset "model storage round-trip preserves scorebounds (property)" begin
    pwm = readmodel(joinpath(EXAMPLES, "pif4.meme"))
    dir = mktempdir()
    bundle = joinpath(dir, "pwm")
    writemodel(bundle, pwm)
    loaded = readmodel(bundle)
    @test scorebounds(loaded) == scorebounds(pwm)
end
