# Tests for the public extension contract (ADR 0003, Extensibility API Plan §4, §9).
#
# Defines a minimal custom model that satisfies only the three required
# methods (`modelname`, `motif_length`, `scan_pair_kernel!`) plus models
# with left/right context, and verifies that the public scan/prepare/
# compare/sites workflows work through the public API only — no private
# names, no struct fields named `representation`/`weights`/`order`/`span`.

using Test
using Mimosa

# ── Minimal model without context ────────────────────────────────────────────
#
# This is a fixed-pattern PWM-like model implemented in pure Julia. It
# produces a deterministic, checkable forward score equal to the number
# of matches to a fixed consensus sequence, and a reverse score equal to
# the number of matches to its reverse complement.

struct ConsensusModel <: Mimosa.AbstractMotifModel
    name::String
    consensus::Vector{UInt8}  # encoded bases (A,C,G,T) — N not allowed in consensus
end

Mimosa.modelname(m::ConsensusModel) = m.name
Mimosa.motif_length(m::ConsensusModel) = length(m.consensus)

function Mimosa.scan_pair_kernel!(
    fwd_out::AbstractVector{Float32},
    rev_out::AbstractVector{Float32},
    model::ConsensusModel,
    sequence::AbstractVector{UInt8},
    n_positions::Int,
)
    consensus = model.consensus
    L = length(consensus)
    rc_consensus = UInt8[b == N_CODE ? b : (0x03 - b) for b in reverse(consensus)]
    @inbounds for pos in 1:n_positions
        fwd = zero(Float32)
        rev = zero(Float32)
        for k in 1:L
            base = sequence[pos + k - 1]
            fwd += (base == consensus[k]) ? 1.0f0 : 0.0f0
            rev += (base == rc_consensus[k]) ? 1.0f0 : 0.0f0
        end
        fwd_out[pos] = fwd
        rev_out[pos] = rev
    end
    return (fwd_out, rev_out)
end

# ── Model with both left and right context ───────────────────────────────────
#
# This model scores the motif site, but its kernel reads `left_context`
# bases before the site and `right_context` bases after it to bias the
# score by the surrounding GC content. The site is still `motif_length`
# bases long.

struct ContextModel <: Mimosa.AbstractMotifModel
    name::String
    motif_length::Int
    left_context::Int
    right_context::Int
end

Mimosa.modelname(m::ContextModel) = m.name
Mimosa.motif_length(m::ContextModel) = m.motif_length
Mimosa.left_context(m::ContextModel) = m.left_context
Mimosa.right_context(m::ContextModel) = m.right_context

function Mimosa.scan_pair_kernel!(
    fwd_out::AbstractVector{Float32},
    rev_out::AbstractVector{Float32},
    model::ContextModel,
    sequence::AbstractVector{UInt8},
    n_positions::Int,
)
    ml = model.motif_length
    lc = model.left_context
    rc = model.right_context
    # The site itself starts at offset `lc` within the window.
    @inbounds for pos in 1:n_positions
        site_start = pos + lc
        # Forward site score = sum of (G or C) inside the motif site.
        fwd = zero(Float32)
        for k in 1:ml
            base = sequence[site_start + k - 1]
            fwd += (base == 0x01 || base == 0x02) ? 1.0f0 : 0.0f0
        end
        # Context contribution: count GC in left+right context.
        gc_ctx = 0.0f0
        for k in 1:lc
            base = sequence[pos + k - 1]
            gc_ctx += (base == 0x01 || base == 0x02) ? 0.5f0 : 0.0f0
        end
        for k in 1:rc
            base = sequence[site_start + ml + k - 1]
            gc_ctx += (base == 0x01 || base == 0x02) ? 0.5f0 : 0.0f0
        end
        fwd_out[pos] = fwd + gc_ctx
        # Reverse orientation: scan the same physical window but reverse-
        # complement the bases. The reverse site is the same physical
        # interval, only the returned bases would be RC'd. We mirror the
        # forward formula on the reverse-complement of the window so the
        # reverse score is well-defined and the kernel is self-consistent.
        rev = zero(Float32)
        window_len = lc + ml + rc
        for k in 1:ml
            j = window_len - (lc + k - 1)  # mirrored site position in RC window
            base = sequence[pos + j - 1]
            base_rc = base == N_CODE ? base : (0x03 - base)
            rev += (base_rc == 0x01 || base_rc == 0x02) ? 1.0f0 : 0.0f0
        end
        rev_out[pos] = rev
    end
    return (fwd_out, rev_out)
end

# ── Helpers ──────────────────────────────────────────────────────────────────

function _encode(seq::AbstractString)
    return Mimosa.encode_sequence(seq)
end

# ── Tests: validate_model ────────────────────────────────────────────────────

@testset "Extension: validate_model for minimal model" begin
    m = ConsensusModel("cm", _encode("ACGT"))
    @test validate_model(m; capability=:compare) === m
    @test validate_model(m; capability=:sites) === m
    # No model_fingerprint defined → :cache must fail with a clear error.
    @test_throws ModelInterfaceError validate_model(m; capability=:cache)
    @test_throws ModelInterfaceError model_collection_fingerprint(AbstractProfileSource[m])
end

@testset "Extension: validate_model for context model" begin
    m = ContextModel("ctx", 4, 2, 1)
    @test validate_model(m; capability=:compare) === m
    @test validate_model(m; capability=:sites) === m
    @test Mimosa.motif_length(m) == 4
    @test Mimosa.left_context(m) == 2
    @test Mimosa.right_context(m) == 1
    @test Mimosa.window_size(m) == 7
    @test Mimosa.site_start_offset(m) == 2
    @test Mimosa.npositions(m, 10) == 4
    @test Mimosa.npositions(m, 7) == 1
    @test Mimosa.npositions(m, 6) == 0
end

@testset "Extension: validate_model rejects bad interface" begin
    # Missing modelname.
    struct NoNameModel <: Mimosa.AbstractMotifModel
        length::Int
    end
    Mimosa.motif_length(m::NoNameModel) = m.length
    Mimosa.scan_pair_kernel!(
        f::AbstractVector{Float32},
        r::AbstractVector{Float32},
        m::NoNameModel,
        s::AbstractVector{UInt8},
        n::Int,
    ) = (f, r)
    @test_throws ModelInterfaceError validate_model(NoNameModel(3); capability=:compare)

    # Missing motif_length.
    struct NoLengthModel <: Mimosa.AbstractMotifModel
        name::String
    end
    Mimosa.modelname(m::NoLengthModel) = m.name
    @test_throws ModelInterfaceError validate_model(NoLengthModel("x"); capability=:compare)

    # Missing scan capability.
    struct NoScanModel <: Mimosa.AbstractMotifModel
        name::String
        length::Int
    end
    Mimosa.modelname(m::NoScanModel) = m.name
    Mimosa.motif_length(m::NoScanModel) = m.length
    @test_throws ModelInterfaceError validate_model(
        NoScanModel("x", 3); capability=:compare
    )

    # Non-positive motif_length.
    struct ZeroLengthModel <: Mimosa.AbstractMotifModel
        name::String
    end
    Mimosa.modelname(m::ZeroLengthModel) = m.name
    Mimosa.motif_length(m::ZeroLengthModel) = 0
    @test_throws ModelInterfaceError validate_model(
        ZeroLengthModel("x"); capability=:compare
    )

    # Empty name.
    struct EmptyNameModel <: Mimosa.AbstractMotifModel
        name::String
        length::Int
    end
    Mimosa.modelname(m::EmptyNameModel) = m.name
    Mimosa.motif_length(m::EmptyNameModel) = m.length
    Mimosa.scan_pair_kernel!(
        f::AbstractVector{Float32},
        r::AbstractVector{Float32},
        m::EmptyNameModel,
        s::AbstractVector{UInt8},
        n::Int,
    ) = (f, r)
    @test_throws ModelInterfaceError validate_model(
        EmptyNameModel("", 3); capability=:compare
    )
end

@testset "Extension: validation boundaries and geometry overrides" begin
    @test !isdefined(Mimosa, :AbstractMatrixMotif)
    @test !isdefined(Mimosa, :AbstractHigherOrderMotif)
    @test_throws ModelInterfaceError validate_model(42; capability=:compare)
    @test_throws ArgumentError validate_model(
        ConsensusModel("cm", _encode("ACGT")); capability=:storage_write
    )

    struct BoundaryNoNameModel <: Mimosa.AbstractMotifModel end
    Mimosa.motif_length(::BoundaryNoNameModel) = 1
    Mimosa.scan_pair_kernel!(
        f::AbstractVector{Float32},
        r::AbstractVector{Float32},
        ::BoundaryNoNameModel,
        ::AbstractVector{UInt8},
        ::Int,
    ) = (f, r)
    @test_throws ModelInterfaceError scan(BoundaryNoNameModel(), UInt8[0])

    struct BadWindowModel <: Mimosa.AbstractMotifModel end
    Mimosa.modelname(::BadWindowModel) = "bad-window"
    Mimosa.motif_length(::BadWindowModel) = 2
    Mimosa.window_size(::BadWindowModel) = 3
    Mimosa.scan_pair_kernel!(
        f::AbstractVector{Float32},
        r::AbstractVector{Float32},
        ::BadWindowModel,
        ::AbstractVector{UInt8},
        ::Int,
    ) = (f, r)
    @test_throws ModelInterfaceError validate_model(BadWindowModel())

    struct BadOffsetModel <: Mimosa.AbstractMotifModel end
    Mimosa.modelname(::BadOffsetModel) = "bad-offset"
    Mimosa.motif_length(::BadOffsetModel) = 2
    Mimosa.site_start_offset(::BadOffsetModel) = 1
    Mimosa.scan_pair_kernel!(
        f::AbstractVector{Float32},
        r::AbstractVector{Float32},
        ::BadOffsetModel,
        ::AbstractVector{UInt8},
        ::Int,
    ) = (f, r)
    @test validate_model(BadOffsetModel(); capability=:compare) isa BadOffsetModel
    @test_throws ModelInterfaceError validate_model(BadOffsetModel(); capability=:sites)
end

@testset "Extension: abstract accessor result types are canonicalized" begin
    struct FlexibleModel{S<:AbstractString,I<:Integer} <: Mimosa.AbstractMotifModel
        name::S
        width::I
    end
    Mimosa.modelname(model::FlexibleModel) = model.name
    Mimosa.motif_length(model::FlexibleModel) = model.width
    function Mimosa.scan_pair_kernel!(
        f::AbstractVector{Float32},
        r::AbstractVector{Float32},
        ::FlexibleModel,
        ::AbstractVector{UInt8},
        n::Int,
    )
        fill!(f, 1.0f0)
        fill!(r, 2.0f0)
        return (f, r)
    end

    substring_model = FlexibleModel(SubString("substring", 1, 3), Int32(1))
    batch = EncodedSequenceBatch([UInt8[0, 1, 2]])
    prepared = prepare_profile(substring_model, batch)
    @test prepared.name == "sub"
    @test prepared.name isa String

    big_model = FlexibleModel("big", big(1))
    @test validate_model(big_model) === big_model
    @test scan(big_model, UInt8[0, 1]) == Float32[1, 1]

    dest = zeros(Float64, 2)
    @test scan!(dest, big_model, UInt8[0, 1]) === dest
    @test dest == [1.0, 1.0]
end

@testset "Extension: pair kernel return and fingerprint validation" begin
    struct BadReturnModel <: Mimosa.AbstractMotifModel end
    Mimosa.modelname(::BadReturnModel) = "bad-return"
    Mimosa.motif_length(::BadReturnModel) = 1
    function Mimosa.scan_pair_kernel!(
        f::AbstractVector{Float32},
        r::AbstractVector{Float32},
        ::BadReturnModel,
        ::AbstractVector{UInt8},
        ::Int,
    )
        fill!(f, 0.0f0)
        fill!(r, 0.0f0)
        return nothing
    end
    @test_throws InvariantError scan(BadReturnModel(), UInt8[0])

    struct WeakFingerprintModel <: Mimosa.AbstractMotifModel end
    Mimosa.modelname(::WeakFingerprintModel) = "weak-fingerprint"
    Mimosa.motif_length(::WeakFingerprintModel) = 1
    Mimosa.scan_pair_kernel!(
        f::AbstractVector{Float32},
        r::AbstractVector{Float32},
        ::WeakFingerprintModel,
        ::AbstractVector{UInt8},
        ::Int,
    ) = (f, r)
    Mimosa.model_fingerprint(::WeakFingerprintModel) = "not-a-sha256"
    @test_throws ModelInterfaceError validate_model(
        WeakFingerprintModel(); capability=:cache
    )

    struct BrokenFingerprintModel <: Mimosa.AbstractMotifModel end
    Mimosa.modelname(::BrokenFingerprintModel) = "broken-fingerprint"
    Mimosa.motif_length(::BrokenFingerprintModel) = 1
    Mimosa.scan_pair_kernel!(
        f::AbstractVector{Float32},
        r::AbstractVector{Float32},
        ::BrokenFingerprintModel,
        ::AbstractVector{UInt8},
        ::Int,
    ) = (f, r)
    Mimosa.model_fingerprint(::BrokenFingerprintModel) = push!(1, 1)
    @test_throws MethodError validate_model(BrokenFingerprintModel(); capability=:cache)
end

# ── Tests: scanning through the public API ────────────────────────────────────

@testset "Extension: scan single sequence" begin
    m = ConsensusModel("cm", _encode("ACGT"))
    seq = _encode("GGGGACGTGGGG")  # one forward ACGT site at position 5
    n_pos = Mimosa.npositions(m, length(seq))
    @test n_pos == length(seq) - 4 + 1

    fwd = scan(m, seq; strands=ForwardOnly())
    rev = scan(m, seq; strands=ReverseOnly())
    best = scan(m, seq; strands=BestStrand())
    both = scan(m, seq; strands=BothStrands())

    @test length(fwd) == n_pos
    @test length(rev) == n_pos
    @test length(best) == n_pos
    @test length(both.forward) == n_pos
    @test length(both.reverse) == n_pos

    # The single ACGT window at position 5 should give score 4 (perfect match).
    @test fwd[5] == 4.0f0
    # No other position should achieve the perfect-match score.
    @test count(==(4.0f0), fwd) == 1

    # Best must equal max of forward/reverse.
    @test best == max.(fwd, rev)

    # BothStrands returns the same forward and reverse tracks.
    @test both.forward == fwd
    @test both.reverse == rev
end

@testset "Extension: scan empty / too-short / exact-window sequences" begin
    m = ConsensusModel("cm", _encode("ACGT"))
    @test scan(m, UInt8[]) == Float32[]
    @test scan(m, _encode("AC")) == Float32[]
    @test length(scan(m, _encode("ACGT"))) == 1
    @test scan(m, _encode("ACGT"); strands=ForwardOnly())[1] == 4.0f0
end

@testset "Extension: scan batch and ragged batches" begin
    m = ConsensusModel("cm", _encode("ACGT"))
    rows = [
        _encode("ACGTACGT"),         # 5 positions
        UInt8[],                     # 0 positions
        _encode("AC"),               # 0 positions (too short)
        _encode("GGGGACGTGGGG"),    # 9 positions, one perfect site
    ]
    batch = EncodedSequenceBatch(rows)

    serial_fwd = scan(m, batch; strands=ForwardOnly(), execution=SerialExecution())
    @test nrows(serial_fwd) == 4
    @test rowlength(serial_fwd, 1) == 5
    @test rowlength(serial_fwd, 2) == 0
    @test rowlength(serial_fwd, 3) == 0
    @test rowlength(serial_fwd, 4) == 9
    @test row(serial_fwd, 4)[5] == 4.0f0

    both_serial = scan(m, batch; strands=BothStrands(), execution=SerialExecution())
    both_threaded = scan(m, batch; strands=BothStrands(), execution=ThreadedExecution(2))
    @test both_serial.forward == both_threaded.forward
    @test both_serial.reverse == both_threaded.reverse

    # BestStrand threaded == serial, order preserved.
    best_serial = scan(m, batch; strands=BestStrand(), execution=SerialExecution())
    best_threaded = scan(m, batch; strands=BestStrand(), execution=ThreadedExecution(2))
    @test best_serial == best_threaded
end

@testset "Extension: context model scan honors geometry" begin
    # motif_length=3, left_context=1, right_context=1 → window_size=5
    m = ContextModel("ctx", 3, 1, 1)
    @test Mimosa.window_size(m) == 5
    seq = _encode("ACGTACGTACGT")  # length 12 → 8 positions
    n_pos = Mimosa.npositions(m, length(seq))
    @test n_pos == 8

    fwd = scan(m, seq; strands=ForwardOnly())
    @test length(fwd) == n_pos
    @test all(isfinite, fwd)
end

# ── Tests: prepare_profile and compare through the public API ────────────────

@testset "Extension: prepare_profile and compare" begin
    m1 = ConsensusModel("query", _encode("ACGT"))
    m2 = ConsensusModel("target", _encode("ACGT"))
    sequences = Mimosa.make_random_sequences(10, 60; seed=42)

    prepared = prepare_profile(m1, sequences)
    @test prepared isa PreparedProfile
    @test Mimosa.modelname(prepared) == "query"

    # Motif-to-motif compare.
    result = compare(m1, m2, sequences; metric=:co, search_range=3, window_radius=2)
    @test result isa ComparisonResult
    @test result.query == "query"
    @test result.target == "target"

    # Prepared query vs motif target.
    res2 = compare(prepared, m2, sequences; metric=:co, search_range=3, window_radius=2)
    @test res2 isa ComparisonResult
    @test res2.query == "query"

    # Motif query vs prepared target.
    prepared2 = prepare_profile(m2, sequences)
    res3 = compare(m1, prepared2, sequences; metric=:co, search_range=3, window_radius=2)
    @test res3 isa ComparisonResult
    @test res3.target == "target"

    # One-to-many.
    targets = AbstractMotifModel[m2, ConsensusModel("other", _encode("ACGA"))]
    results = compare(
        prepared, targets, sequences; metric=:co, search_range=3, window_radius=2
    )
    @test length(results) == 2
    @test all(r isa ComparisonResult for r in results)
end

@testset "Extension: custom/built-in comparison in both orders" begin
    custom = ConsensusModel("custom", _encode("ACGT"))
    pwm = readmodel(joinpath(dirname(dirname(@__DIR__)), "examples", "pif4.meme"))
    sequences = Mimosa.make_random_sequences(8, 80; seed=11)

    r1 = compare(custom, pwm, sequences; metric=:co, search_range=4, window_radius=2)
    r2 = compare(pwm, custom, sequences; metric=:co, search_range=4, window_radius=2)
    @test r1.query == "custom"
    @test r1.target == pwm.name
    @test r2.query == pwm.name
    @test r2.target == "custom"
end

# ── Tests: sites and PFM reconstruction through the public API ───────────────

@testset "Extension: selectsites and reconstruct_pfm" begin
    m = ConsensusModel("cm", _encode("ACGT"))
    batch = Mimosa.make_random_sequences(20, 100; seed=7)

    sites = selectsites(m, batch, BestPerSequence(); strands=BothStrands())
    @test sites isa SiteCollection

    # The start of each site must be the scan position. The reconstructed
    # site window length must equal motif_length.
    if !isempty(sites)
        pfm = reconstruct_pfm(m, batch, BestPerSequence(); pseudocount=Float32(1e-4))
        @test size(pfm) == (4, 4)
        @test all(isfinite, pfm)
        # Each column of the PFM must sum to ~1 (with pseudocount).
        @test all(isapprox.(sum(pfm; dims=1), 1.0f0; atol=1e-3))
    end
end

@testset "Extension: context model sites extract only the site window" begin
    # motif_length=3, left_context=1 → site_start_offset=1
    m = ContextModel("ctx", 3, 1, 0)
    @test Mimosa.site_start_offset(m) == 1
    # Build a batch where we know the site boundaries.
    batch = EncodedSequenceBatch([_encode("GACGTAAAA")])
    sites = selectsites(m, batch, BestPerSequence(); strands=BothStrands())
    @test !isempty(sites)
    pfm = reconstruct_pfm(m, batch, BestPerSequence(); pseudocount=0.25f0)
    @test size(pfm) == (4, 3)
end

# ── Tests: cache capability without fingerprint ──────────────────────────────

@testset "Extension: compare does not require fingerprint" begin
    m1 = ConsensusModel("q", _encode("ACGT"))
    m2 = ConsensusModel("t", _encode("ACGT"))
    sequences = Mimosa.make_random_sequences(5, 40; seed=3)
    # No fingerprint method defined on ConsensusModel; compare must still work.
    @test_nowarn compare(m1, m2, sequences; metric=:co, search_range=2, window_radius=1)
end

@testset "Extension: fingerprint capability diagnostic" begin
    m = ConsensusModel("cm", _encode("ACGT"))
    # validate_model(:cache) must fail with a ModelInterfaceError that names
    # the fingerprint method.
    err = try
        validate_model(m; capability=:cache)
        @test false  # should have thrown
    catch e
        e
    end
    @test err isa ModelInterfaceError
    @test occursin("model_fingerprint", err.message)
end

# ── Tests: worker exception propagation ──────────────────────────────────────

@testset "Extension: worker kernel exception propagates" begin
    struct ThrowingModel <: Mimosa.AbstractMotifModel
        name::String
        length::Int
    end
    Mimosa.modelname(m::ThrowingModel) = m.name
    Mimosa.motif_length(m::ThrowingModel) = m.length
    function Mimosa.scan_pair_kernel!(
        f::AbstractVector{Float32},
        r::AbstractVector{Float32},
        m::ThrowingModel,
        s::AbstractVector{UInt8},
        n::Int,
    )
        return error("kernel failure for $(m.name)")
    end

    m = ThrowingModel("thrower", 4)
    batch = EncodedSequenceBatch([_encode("ACGTACGT"), _encode("GGGGCCCC")])
    err = try
        scan(m, batch; strands=ForwardOnly())
        nothing
    catch e
        e
    end
    @test err !== nothing
    @test occursin("kernel failure", sprint(showerror, err))
end
