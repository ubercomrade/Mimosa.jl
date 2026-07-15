using Test
using Mimosa

const REPO_ROOT = dirname(dirname(@__DIR__))
const EXAMPLES = joinpath(REPO_ROOT, "examples")

@testset "encoding" begin
    @test encode_sequence("ACGT") == UInt8[0x00, 0x01, 0x02, 0x03]
    @test encode_sequence("acgt") == UInt8[0x00, 0x01, 0x02, 0x03]
    @test encode_sequence("ACGTACGT") ==
        UInt8[0x00, 0x01, 0x02, 0x03, 0x00, 0x01, 0x02, 0x03]
    @test encode_sequence("N") == UInt8[0x04]
    @test encode_sequence("NNNN") == UInt8[0x04, 0x04, 0x04, 0x04]
    @test encode_sequence("") == UInt8[]
    # IUPAC and other non-ACGT map to N
    @test encode_sequence("RYWSKMBDHV") == UInt8[fill(0x04, 10)...]
    # Mixed case and special chars
    @test encode_sequence("AcGt") == UInt8[0x00, 0x01, 0x02, 0x03]
    @test encode_sequence("AC-TG") == UInt8[0x00, 0x01, 0x04, 0x03, 0x02]

    # encode_base
    @test encode_base(UInt8('A')) == 0x00
    @test encode_base(UInt8('T')) == 0x03
    @test encode_base(UInt8('N')) == 0x04
    @test encode_base(UInt8('-')) == 0x04
end

@testset "Sequence layout and model capability contracts" begin
    pwm = PWM("p", ones(Float32, 5, 2), (0.25f0, 0.25f0, 0.25f0, 0.25f0))
    @test is_scannable(pwm)
    @test_throws ArgumentError from_padded(UInt8[0 1; 2 3], Int[1])
    float_pwm = PWM("float", ones(Float64, 5, 2), (0.25, 0.25, 0.25, 0.25))
    @test scan(float_pwm, UInt8[0, 1, 2]) isa Vector{Float32}
end

@testset "reverse_complement" begin
    # Basic complement
    @test reverse_complement(UInt8[0x00]) == UInt8[0x03]  # A -> T
    @test reverse_complement(UInt8[0x01]) == UInt8[0x02]  # C -> G
    @test reverse_complement(UInt8[0x02]) == UInt8[0x01]  # G -> C
    @test reverse_complement(UInt8[0x03]) == UInt8[0x00]  # T -> A
    @test reverse_complement(UInt8[0x04]) == UInt8[0x04]  # N -> N

    # Reverse: ACGT -> ACGT (palindrome in 5-ary: complement of ACGT = TGC A, reversed = ACGT)
    @test reverse_complement(encode_sequence("ACGT")) == encode_sequence("ACGT")
    @test reverse_complement(encode_sequence("AAAA")) == encode_sequence("TTTT")
    @test reverse_complement(encode_sequence("ACGTACGT")) == encode_sequence("ACGTACGT")

    # Empty
    @test reverse_complement(UInt8[]) == UInt8[]

    # In-place variant: ACGT is a reverse-complement palindrome
    src = UInt8[0x00, 0x01, 0x02, 0x03]
    dest = similar(src)
    reverse_complement!(dest, src)
    @test dest == UInt8[0x00, 0x01, 0x02, 0x03]  # ACGT is RC palindrome
    # src unchanged
    @test src == UInt8[0x00, 0x01, 0x02, 0x03]
    # Non-palindrome: AAAA → TTTT
    src2 = UInt8[0x00, 0x00, 0x00, 0x00]
    dest2 = similar(src2)
    reverse_complement!(dest2, src2)
    @test dest2 == UInt8[0x03, 0x03, 0x03, 0x03]
end

@testset "EncodedSequenceBatch" begin
    # Basic construction
    rows = [UInt8[0, 1, 2, 3], UInt8[0, 1, 2, 3, 0, 1, 2, 3], UInt8[]]
    batch = EncodedSequenceBatch(rows)
    @test nsequences(batch) == 3
    @test seqlength(batch, 1) == 4
    @test seqlength(batch, 2) == 8
    @test seqlength(batch, 3) == 0
    @test sequence(batch, 1) == UInt8[0, 1, 2, 3]
    @test sequence(batch, 3) == UInt8[]

    # Empty batch
    empty_batch = empty_sequence_batch()
    @test nsequences(empty_batch) == 0

    # Constructor invariants
    @test_throws ArgumentError EncodedSequenceBatch(UInt8[0, 1], Int[0, 3])  # offsets[1] != 1
    @test_throws ArgumentError EncodedSequenceBatch(UInt8[0, 1], Int[1, 0, 3])  # not monotonic

    # Iteration
    batch2 = EncodedSequenceBatch([UInt8[0, 1], UInt8[2, 3]])
    seqs = collect(batch2)
    @test length(seqs) == 2
    @test seqs[1] == UInt8[0, 1]
    @test seqs[2] == UInt8[2, 3]

    # to_padded / from_padded round-trip
    padded, lengths = to_padded(batch)
    @test size(padded) == (3, 8)
    @test lengths == [4, 8, 0]
    @test padded[1, 1:4] == UInt8[0, 1, 2, 3]
    @test padded[3, 1] == N_CODE  # empty seq padded
    rt = from_padded(padded, lengths)
    @test rt == batch
end

@testset "RaggedArray" begin
    # Basic
    data = Float32[1.0, 2.0, 3.0, 4.0, 5.0]
    offsets = Int[1, 3, 5, 6]  # rows: [1,2], [3,4], [5]
    rag = RaggedArray(data, offsets)
    @test nrows(rag) == 3
    @test rowlength(rag, 1) == 2
    @test rowlength(rag, 2) == 2
    @test rowlength(rag, 3) == 1
    @test row(rag, 1) == Float32[1.0, 2.0]
    @test row(rag, 3) == Float32[5.0]

    # build_ragged
    rag2 = build_ragged([[1.0f0, 2.0f0], [3.0f0], Float32[]])
    @test nrows(rag2) == 3
    @test rowlength(rag2, 3) == 0

    # Empty
    empty_rag = empty_ragged(Float32)
    @test nrows(empty_rag) == 0

    # Validation
    @test_throws ArgumentError RaggedArray(Float32[1.0], Int[2])  # offsets[1] != 1
    @test_throws ArgumentError RaggedArray(Float32[1.0], Int[1, 3])  # offsets[end] != len+1
end

@testset "FASTA reader" begin
    batch, names = read_fasta(joinpath(EXAMPLES, "foreground.fa"))
    @test nsequences(batch) == 100
    @test startswith(names[1], "peaks_0")
    # All lengths positive
    for i in 1:nsequences(batch)
        @test seqlength(batch, i) > 0
    end

    # Padded conversion matches Python format
    padded, lengths = to_padded(batch)
    @test size(padded, 1) == 100
    @test all(lengths .> 0)

    # Empty FASTA
    mktemp() do path, io
        write(io, "")
        close(io)
        @test_throws ModelFormatError read_fasta(path)
    end

    # FASTA with empty sequence
    mktemp() do path, io
        write(io, ">seq1\nACGT\n>empty\n>seq2\nGGGG\n")
        close(io)
        b, n = read_fasta(path)
        @test nsequences(b) == 3
        @test seqlength(b, 2) == 0  # empty sequence
        @test n == ["seq1", "empty", "seq2"]
    end
end

@testset "PWM scan single sequence" begin
    pwm = readmodel(joinpath(EXAMPLES, "pif4.meme"))
    W = length(pwm)
    @test W == 12

    # Short sequence (shorter than motif)
    short_seq = encode_sequence("ACGT")
    @test scan(pwm, short_seq; strands=ForwardOnly()) == Float32[]
    @test scan(pwm, short_seq; strands=ReverseOnly()) == Float32[]
    @test scan(pwm, short_seq; strands=BestStrand()) == Float32[]
    pair = scan(pwm, short_seq; strands=BothStrands())
    @test pair.forward == Float32[]
    @test pair.reverse == Float32[]

    # All-N sequence
    n_seq = encode_sequence("NNNNNNNNNNNNNN")
    fwd = scan(pwm, n_seq; strands=ForwardOnly())
    @test length(fwd) == 3  # 14 - 12 + 1
    # N positions should use the minimum row (row 5)
    @test all(isfinite, fwd)

    # Normal sequence
    seq = encode_sequence("ACGTACGTACGTACGTACGTACGTAC")
    fwd = scan(pwm, seq; strands=ForwardOnly())
    rev = scan(pwm, seq; strands=ReverseOnly())
    best = scan(pwm, seq; strands=BestStrand())
    @test length(fwd) == length(seq) - W + 1
    @test length(rev) == length(fwd)
    @test length(best) == length(fwd)

    # Best strand is element-wise max of forward and reverse
    @test best == max.(fwd, rev)

    # Both strands
    pair = scan(pwm, seq; strands=BothStrands())
    @test pair.forward == fwd
    @test pair.reverse == rev
end

@testset "PWM scan in-place" begin
    pwm = readmodel(joinpath(EXAMPLES, "pif4.meme"))
    seq = encode_sequence("ACGTACGTACGTACGTACGTACGTAC")

    # Forward
    fwd_alloc = scan(pwm, seq; strands=ForwardOnly())
    dest = Vector{Float32}(undef, length(fwd_alloc))
    scan!(dest, pwm, seq; strands=ForwardOnly())
    @test dest == fwd_alloc

    # Reverse
    rev_alloc = scan(pwm, seq; strands=ReverseOnly())
    dest2 = Vector{Float32}(undef, length(rev_alloc))
    scan!(dest2, pwm, seq; strands=ReverseOnly())
    @test dest2 == rev_alloc

    # Best
    best_alloc = scan(pwm, seq; strands=BestStrand())
    dest3 = Vector{Float32}(undef, length(best_alloc))
    scan!(dest3, pwm, seq; strands=BestStrand())
    @test dest3 == best_alloc

    # Destination too small
    @test_throws ArgumentError scan!(Vector{Float32}(undef, 1), pwm, seq)

    # BothStrands with scan! is not supported
    @test_throws ArgumentError scan!(
        Vector{Float32}(undef, 10), pwm, seq; strands=BothStrands()
    )
end

@testset "PWM scan batch" begin
    pwm = readmodel(joinpath(EXAMPLES, "pif4.meme"))
    W = length(pwm)

    batch = EncodedSequenceBatch([
        encode_sequence("ACGTACGTACGTACGTACGTACGTAC"),
        encode_sequence("TTTTGGGGCCCCAAAATTTTGGGGCCC"),
        encode_sequence("ACGT"),  # short
    ])

    # Forward
    fwd = scan(pwm, batch; strands=ForwardOnly())
    @test nrows(fwd) == 3
    seq1_len = length(encode_sequence("ACGTACGTACGTACGTACGTACGTAC"))
    @test rowlength(fwd, 1) == seq1_len - W + 1
    @test rowlength(fwd, 3) == 0  # short sequence

    # Verify batch == single-sequence
    for i in 1:nsequences(batch)
        single = scan(pwm, sequence(batch, i); strands=ForwardOnly())
        @test row(fwd, i) == single
    end

    # Reverse
    rev = scan(pwm, batch; strands=ReverseOnly())
    @test nrows(rev) == 3

    # Best
    best = scan(pwm, batch; strands=BestStrand())
    for i in 1:nsequences(batch)
        @test row(best, i) == max.(row(fwd, i), row(rev, i))
    end

    # Both strands
    pair = scan(pwm, batch; strands=BothStrands())
    @test pair.forward == fwd
    @test pair.reverse == rev
end

@testset "scan score bounds" begin
    pwm = readmodel(joinpath(EXAMPLES, "pif4.meme"))
    mn, mx = scorebounds(pwm)

    # No scan score should exceed theoretical bounds
    seq = encode_sequence("ACGTACGTACGTACGTACGTACGTAC")
    fwd = scan(pwm, seq; strands=ForwardOnly())
    @test all(s -> mn - 1e-3 <= s <= mx + 1e-3, fwd)
    rev = scan(pwm, seq; strands=ReverseOnly())
    @test all(s -> mn - 1e-3 <= s <= mx + 1e-3, rev)
end

@testset "StrandPair" begin
    pair = StrandPair([1.0f0, 2.0f0], [3.0f0, 4.0f0])
    @test pair.forward == [1.0f0, 2.0f0]
    @test pair.reverse == [3.0f0, 4.0f0]
end
