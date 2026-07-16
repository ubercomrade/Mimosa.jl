# Downstream Contract Test for Mimosa.jl
#
# This test verifies that a downstream consumer (e.g., MotifHORDE.jl) can
# use Mimosa.jl through only the documented public API, without importing
# internal submodules or accessing internal functions.
#
# Run with:
#   julia --project=test/downstream test/downstream/runtests.jl
#
# The test environment is independent of the main test suite: it does not
# include any shared test helpers, NPY readers, or fixture metadata.

using Pkg
# Ensure Mimosa is available from the local package root.
Pkg.develop(; path=dirname(dirname(@__DIR__)))

using Test
using Mimosa

# Path to example data files (repo-level, not test internals)
const REPO_ROOT = dirname(dirname(@__DIR__))
const EXAMPLES = joinpath(REPO_ROOT, "examples")

# ---------------------------------------------------------------------------
# Verify that all exported names are accessible
# ---------------------------------------------------------------------------
@testset "Downstream contract: exports" begin
    # Model I/O
    @test isdefined(Mimosa, :readmodel)
    @test isdefined(Mimosa, :writemodel)
    @test isdefined(Mimosa, :readsequences)

    # Scanning
    @test isdefined(Mimosa, :scan)
    @test isdefined(Mimosa, :scan!)
    @test isdefined(Mimosa, :scorebounds)

    # Comparison
    @test isdefined(Mimosa, :compare)
    @test isdefined(Mimosa, :ComparisonResult)

    # Site extraction
    @test isdefined(Mimosa, :selectsites)
    @test isdefined(Mimosa, :reconstruct_pfm)

    # Statistics
    @test isdefined(Mimosa, :build_null)
    @test isdefined(Mimosa, :pvalue)
    @test isdefined(Mimosa, :adjusted_pvalues)
    @test isdefined(Mimosa, :evalue)
    @test isdefined(Mimosa, :savenull)
    @test isdefined(Mimosa, :loadnull)

    # Execution policies
    @test isdefined(Mimosa, :SerialExecution)
    @test isdefined(Mimosa, :ThreadedExecution)

    # Cache
    @test isdefined(Mimosa, :Cache)
    @test isdefined(Mimosa, :clearcache)

    # Serialization
    @test isdefined(Mimosa, :to_json)
    @test isdefined(Mimosa, :to_dict)

    # Errors
    @test isdefined(Mimosa, :MimosaError)
    @test isdefined(Mimosa, :ModelFormatError)
    @test isdefined(Mimosa, :ModelDimensionError)
    @test isdefined(Mimosa, :InvariantError)

    # Annotation
    @test isdefined(Mimosa, :annotate_results)
    @test isdefined(Mimosa, :AnnotatedResult)
    @test isdefined(Mimosa, :ANNOTATED_RESULT_SCHEMA_VERSION)

    # Null distribution types
    @test isdefined(Mimosa, :NullDistribution)
    @test isdefined(Mimosa, :NullBuildConfig)
    @test isdefined(Mimosa, :NullBuildResult)

    # Storage format versions
    @test isdefined(Mimosa, :MODEL_FORMAT_VERSION)
    @test isdefined(Mimosa, :NULL_FORMAT_VERSION)
end

# ---------------------------------------------------------------------------
# Model I/O: read models from files and write/read bundles
# ---------------------------------------------------------------------------
@testset "Downstream contract: model I/O" begin
    # Read PWM from MEME format (readmodel converts PFM→PWM internally)
    pwm = readmodel(joinpath(EXAMPLES, "pif4.meme"))
    @test pwm isa PWM

    # Read PFM file (returns PWM after conversion)
    pwm_from_pfm = readmodel(joinpath(EXAMPLES, "pif4.pfm"))
    @test pwm_from_pfm isa PWM

    # Read BaMM
    bamm = readmodel(joinpath(EXAMPLES, "foxa2.ihbcp"))
    @test bamm isa BaMM

    # Write and re-read a portable bundle (signature: writemodel(path, model))
    tmpdir = mktempdir()
    bundle_path = joinpath(tmpdir, "pwm_bundle")
    writemodel(bundle_path, pwm)
    @test isdir(bundle_path)
    @test isfile(joinpath(bundle_path, "manifest.toml"))

    loaded = readmodel(bundle_path)
    @test loaded isa PWM
    @test loaded.name == pwm.name
    @test size(loaded.representation) == size(pwm.representation)

    # Write and re-read BaMM bundle
    bamm_path = joinpath(tmpdir, "bamm_bundle")
    writemodel(bamm_path, bamm)
    loaded_bamm = readmodel(bamm_path)
    @test loaded_bamm isa BaMM
    @test loaded_bamm.order == bamm.order
    @test loaded_bamm.motif_length == bamm.motif_length
end

# ---------------------------------------------------------------------------
# Sequence reading
# ---------------------------------------------------------------------------
@testset "Downstream contract: sequence reading" begin
    batch, names = readsequences(joinpath(EXAMPLES, "foreground.fa"))
    @test batch isa EncodedSequenceBatch
    @test nsequences(batch) > 0
    @test length(names) == nsequences(batch)
    @test all(seqlength(batch, i) > 0 for i in 1:nsequences(batch))
end

# ---------------------------------------------------------------------------
# Scanning
# ---------------------------------------------------------------------------
@testset "Downstream contract: scanning" begin
    pwm = readmodel(joinpath(EXAMPLES, "pif4.meme"))
    batch, _names = readsequences(joinpath(EXAMPLES, "foreground.fa"))

    # Score bounds
    lo, hi = scorebounds(pwm)
    @test lo isa Float32
    @test hi isa Float32
    @test lo <= hi

    # Scan with different strand policies
    scores_fwd = scan(pwm, batch; strands=ForwardOnly())
    @test scores_fwd isa RaggedArray

    scores_best = scan(pwm, batch; strands=BestStrand())
    @test scores_best isa RaggedArray

    # BothStrands returns a StrandPair{RaggedArray}
    scores_both = scan(pwm, batch; strands=BothStrands())
    @test scores_both isa StrandPair
    @test scores_both.forward isa RaggedArray
    @test scores_both.reverse isa RaggedArray

    # Threaded scan == serial scan
    scores_thr = scan(pwm, batch; strands=BestStrand(), execution=ThreadedExecution(2))
    @test scores_best == scores_thr

    # In-place scan on a single sequence
    seq = sequence(batch, 1)
    n_pos = npositions(seqlength(batch, 1), length(pwm))
    dest = Vector{Float32}(undef, n_pos)
    scan!(dest, pwm, seq; strands=ForwardOnly())
    @test length(dest) == n_pos
    @test all(isfinite, dest)
end

# ---------------------------------------------------------------------------
# Direct motif comparison
# ---------------------------------------------------------------------------
@testset "Downstream contract: profile comparison" begin
    pwm1 = readmodel(joinpath(EXAMPLES, "pif4.meme"))
    pwm2 = readmodel(joinpath(EXAMPLES, "foxa2.meme"))
    sequences = make_random_sequences(4, 80; seed=12)

    last_result = nothing
    for metric in (:co, :dice, :cosine)
        result = compare(pwm1, pwm2, sequences; metric=metric)
        @test result isa ComparisonResult
        @test result.query == pwm1.name
        @test result.target == pwm2.name
        @test result.metric == string(metric)
        last_result = result
    end

    # Self-comparison should give a high score
    self_result = compare(pwm1, pwm1, sequences; metric=:co)
    @test self_result.score >= last_result.score

    # Serialization
    json_str = to_json(last_result)
    @test json_str isa String
    dict = to_dict(last_result)
    @test dict isa Dict
end

# ---------------------------------------------------------------------------
# Profile comparison: one-to-one and one-to-many
# ---------------------------------------------------------------------------
@testset "Downstream contract: profile comparison" begin
    pwm1 = readmodel(joinpath(EXAMPLES, "pif4.meme"))
    pwm2 = readmodel(joinpath(EXAMPLES, "pif4.meme"))

    batch = make_random_sequences(20, 100; seed=42)

    # Motif-derived profile comparison
    result = compare(pwm1, pwm2, batch; metric=:co, search_range=5, window_radius=5)
    @test result isa ComparisonResult
    @test result.n_sites >= 0

    # One-to-many via prepared profile
    sp1 = ScoreProfile("q", scan(pwm1, batch; strands=BestStrand()))
    sp2 = ScoreProfile("t1", scan(pwm2, batch; strands=BestStrand()))
    sp3 = ScoreProfile("t2", scan(pwm1, batch; strands=BestStrand()))

    prepared = prepare_profile(sp1)
    @test prepared isa PreparedProfile

    results = compare(prepared, [sp2, sp3]; metric=:co, search_range=3, window_radius=2)
    @test length(results) == 2
    @test all(r isa ComparisonResult for r in results)
end

# ---------------------------------------------------------------------------
# Site extraction and PFM reconstruction
# ---------------------------------------------------------------------------
@testset "Downstream contract: sites and PFM reconstruction" begin
    pwm = readmodel(joinpath(EXAMPLES, "pif4.meme"))
    batch = make_random_sequences(30, 200; seed=42)

    sites = selectsites(pwm, batch, BestPerSequence(); strands=BothStrands())
    @test sites isa SiteCollection

    pfm = reconstruct_pfm(pwm, batch, BestPerSequence(); pseudocount=Float32(1e-4))
    @test pfm isa AbstractMatrix{Float32}
    @test size(pfm, 1) == 4  # A, C, G, T rows

    # Higher-order model site extraction
    bamm = readmodel(joinpath(EXAMPLES, "foxa2.ihbcp"))
    @test site_start_offset(bamm) == bamm.order

    ho_sites = selectsites(bamm, batch, BestPerSequence(); strands=BothStrands())
    @test ho_sites isa SiteCollection

    if length(ho_sites) > 0
        ho_pfm = reconstruct_pfm(bamm, batch, BestPerSequence(); pseudocount=0.1f0)
        @test ho_pfm isa AbstractMatrix{Float32}
        @test size(ho_pfm, 1) == 4
    end
end

# ---------------------------------------------------------------------------
# Null distribution: build, save, load
# ---------------------------------------------------------------------------
@testset "Downstream contract: null distributions" begin
    pwm1 = readmodel(joinpath(EXAMPLES, "pif4.meme"))
    pwm2 = readmodel(joinpath(EXAMPLES, "foxa2.meme"))

    relations_str = "motif\tgroup\n$(pwm1.name)\tA\n$(pwm2.name)\tB\n"
    rel_path = joinpath(mktempdir(), "groups.tsv")
    write(rel_path, relations_str)
    relations = parse_group_relations(rel_path)

    models = [pwm1, pwm2]
    sequences = make_random_sequences(4, 80; seed=13)
    result = build_null(models, relations; sequences=sequences, metric=:co)
    @test result isa NullBuildResult
    dist = result.distribution
    @test dist isa NullDistribution
    @test dist.strategy == "profile"
    @test dist.metric == "co"
    @test dist.model_collection_fingerprint !== nothing
    @test dist.relation_fingerprint !== nothing
    @test dist.sequence_fingerprint == sequence_fingerprint(sequences)
    @test dist.background_fingerprint == "none"

    # Save and reload null (signature: savenull(path, dist))
    tmpdir = mktempdir()
    null_path = joinpath(tmpdir, "null_dist")
    savenull(null_path, dist)
    @test isdir(null_path)
    @test isfile(joinpath(null_path, "manifest.toml"))

    loaded_dist = loadnull(null_path)
    @test loaded_dist isa NullDistribution
    @test loaded_dist.strategy == dist.strategy
    @test loaded_dist.metric == dist.metric
    @test loaded_dist.n_null == dist.n_null
    @test loaded_dist.model_collection_fingerprint == dist.model_collection_fingerprint
end

# ---------------------------------------------------------------------------
# Annotation
# ---------------------------------------------------------------------------
@testset "Downstream contract: annotation" begin
    # Create multiple PWM models for a meaningful null distribution
    # (GEV fit requires at least 3 scores)
    weights_a = Float32[
        0.5 -0.3 0.8 -0.2 0.1 0.6
        -0.2 0.7 -0.5 0.3 0.8 -0.1
        0.1 -0.4 0.2 0.6 -0.3 0.5
        0.3 0.1 -0.6 0.4 0.2 -0.5
        -0.2 -0.3 -0.5 -0.2 -0.3 -0.5
    ]
    weights_b = Float32[
        0.3 0.2 -0.1 0.5 -0.4
        0.1 0.6 0.3 -0.2 0.7
        -0.4 0.1 0.4 0.3 -0.5
        0.2 -0.3 0.1 -0.1 0.8
        -0.3 -0.4 -0.2 -0.5 -0.6
    ]
    weights_c = Float32[
        -0.1 0.4 0.3 -0.5 0.2
        0.5 -0.2 0.6 0.1 -0.3
        0.3 0.5 -0.4 -0.2 0.7
        -0.5 0.1 0.2 0.8 -0.4
        -0.4 -0.3 -0.5 -0.1 -0.6
    ]
    weights_d = Float32[
        0.4 -0.1 0.7 -0.3 0.5
        -0.3 0.2 0.1 0.6 -0.4
        0.2 -0.5 0.4 -0.1 0.3
        -0.1 0.8 -0.2 -0.4 0.1
        -0.5 -0.4 -0.3 -0.6 -0.2
    ]
    bg = (Float32(0.25), Float32(0.25), Float32(0.25), Float32(0.25))
    models = [
        PWM("m1", weights_a, bg),
        PWM("m2", weights_b, bg),
        PWM("m3", weights_c, bg),
        PWM("m4", weights_d, bg),
    ]

    relations_str = "motif\tgroup\nm1\tA\nm2\tB\nm3\tA\nm4\tB\n"
    rel_path = joinpath(mktempdir(), "groups.tsv")
    write(rel_path, relations_str)
    relations = parse_group_relations(rel_path)

    sequences = make_random_sequences(6, 80; seed=14)
    null_result = build_null(models, relations; sequences=sequences, metric=:co)
    dist = null_result.distribution
    @test dist.n_null >= 3  # GEV fit requires at least 3 scores

    # Create comparison results to annotate
    results = [
        compare(models[1], models[2], sequences; metric=:co),
        compare(models[1], models[3], sequences; metric=:co),
        compare(models[2], models[4], sequences; metric=:co),
    ]

    annotated = annotate_results(results, dist)
    @test length(annotated) == 3
    @test all(a isa AnnotatedResult for a in annotated)
    @test annotated[1].p_value !== nothing
    @test annotated[1].adj_p_value !== nothing
    @test annotated[1].e_value !== nothing
    @test annotated[1].null_n == dist.n_null
    @test annotated[1].null_estimator == "genextreme"
    @test annotated[1].null_id !== nothing

    # With effective_number_of_targets
    annotated2 = annotate_results(results, dist; effective_number_of_targets=50)
    @test annotated2[1].e_value ≈ annotated[1].p_value * 50

    # Serialization
    json_str = to_json(annotated[1])
    @test json_str isa String
    dict = to_dict(annotated[1])
    @test dict isa Dict
end

# ---------------------------------------------------------------------------
# Serial vs threaded equivalence
# ---------------------------------------------------------------------------
@testset "Downstream contract: serial == threaded" begin
    pwm = readmodel(joinpath(EXAMPLES, "pif4.meme"))
    batch = make_random_sequences(20, 200; seed=123)

    serial = scan(pwm, batch; strands=BestStrand(), execution=SerialExecution())
    threaded = scan(pwm, batch; strands=BestStrand(), execution=ThreadedExecution(4))
    @test serial == threaded
end

# ---------------------------------------------------------------------------
# No internal access needed
# ---------------------------------------------------------------------------
@testset "Downstream contract: no internal access needed" begin
    # The fact that all above tests pass with only `using Mimosa`
    # (not `using Mimosa.Scanning` etc.) is the contract.
    @test true
end

# ---------------------------------------------------------------------------
# Extensibility API (ADR 0003): custom model through the public API only
# ---------------------------------------------------------------------------
#
# This testset defines a downstream custom model in a separate module that
# imports Mimosa as a regular dependency. It deliberately avoids any field
# named `representation`, `weights`, `order`, or `span`, and never references
# a `Mimosa._private_name`. The model implements only the three required
# methods plus a left_context override.

module DownstreamCustomModel

using Test
using Mimosa

# A minimal custom model. Note the absence of any field named
# `representation`, `weights`, `order`, or `span`.
struct MatchCounter <: Mimosa.AbstractMotifModel
    label::String          # not named `name`
    pattern::Vector{UInt8} # encoded consensus
end

Mimosa.modelname(m::MatchCounter) = m.label
Mimosa.motif_length(m::MatchCounter) = length(m.pattern)

function Mimosa.scan_kernel!(
    fwd_out::AbstractVector{Float32},
    rev_out::AbstractVector{Float32},
    model::MatchCounter,
    seq::AbstractVector{UInt8},
    n_positions::Int,
)
    pat = model.pattern
    L = length(pat)
    rc_pat = UInt8[b == 0x04 ? b : (0x03 - b) for b in reverse(pat)]
    @inbounds for pos in 1:n_positions
        f = zero(Float32)
        r = zero(Float32)
        for k in 1:L
            b = seq[pos + k - 1]
            f += (b == pat[k]) ? 1.0f0 : 0.0f0
            r += (b == rc_pat[k]) ? 1.0f0 : 0.0f0
        end
        fwd_out[pos] = f
        rev_out[pos] = r
    end
    return (fwd_out, rev_out)
end

struct ContextCounter <: Mimosa.AbstractMotifModel
    label::String
    width::Int
    upstream::Int
    downstream::Int
end

Mimosa.modelname(model::ContextCounter) = model.label
Mimosa.motif_length(model::ContextCounter) = model.width
Mimosa.left_context(model::ContextCounter) = model.upstream
Mimosa.right_context(model::ContextCounter) = model.downstream

function Mimosa.scan_kernel!(
    fwd_out::AbstractVector{Float32},
    rev_out::AbstractVector{Float32},
    model::ContextCounter,
    seq::AbstractVector{UInt8},
    n_positions::Int,
)
    @inbounds for pos in 1:n_positions
        site_start = pos + model.upstream
        score = zero(Float32)
        for offset in 0:(model.width - 1)
            score += seq[site_start + offset] < 0x04 ? 1.0f0 : 0.0f0
        end
        fwd_out[pos] = score
        rev_out[pos] = score
    end
    return (fwd_out, rev_out)
end

function run()
    @testset "validate_model" begin
        m = MatchCounter("downstream", Mimosa.encode_sequence("ACGT"))
        @test validate_model(m; capability=:compare) === m
        @test validate_model(m; capability=:sites) === m
        @test_throws ModelInterfaceError validate_model(m; capability=:cache)
    end

    @testset "scan, batch, and threaded equivalence" begin
        m = MatchCounter("downstream", Mimosa.encode_sequence("ACGT"))
        seq = Mimosa.encode_sequence("GGGGACGTGGGGACGTGGGG")
        n_pos = Mimosa.npositions(m, length(seq))
        @test n_pos == length(seq) - 4 + 1

        fwd = scan(m, seq; strands=ForwardOnly())
        rev = scan(m, seq; strands=ReverseOnly())
        best = scan(m, seq; strands=BestStrand())
        both = scan(m, seq; strands=BothStrands())
        @test length(fwd) == n_pos
        @test best == max.(fwd, rev)
        @test both.forward == fwd
        @test both.reverse == rev

        rows = [
            Mimosa.encode_sequence("ACGTACGT"),
            UInt8[],
            Mimosa.encode_sequence("GGGGACGTGGGG"),
        ]
        batch = EncodedSequenceBatch(rows)
        s_fwd = scan(m, batch; strands=ForwardOnly(), execution=SerialExecution())
        t_fwd = scan(m, batch; strands=ForwardOnly(), execution=ThreadedExecution(2))
        @test s_fwd == t_fwd
        @test nrows(s_fwd) == 3
        @test rowlength(s_fwd, 2) == 0
    end

    @testset "asymmetric context geometry and sites" begin
        model = ContextCounter("context", 3, 1, 2)
        @test validate_model(model; capability=:sites) === model
        @test window_size(model) == 6
        @test site_start_offset(model) == 1

        batch = EncodedSequenceBatch([Mimosa.encode_sequence("AACGTACGTA")])
        scores = scan(model, batch; strands=BothStrands())
        @test rowlength(scores.forward, 1) == 5
        sites = selectsites(model, batch, BestPerSequence())
        @test !isempty(sites)
        pfm = reconstruct_pfm(model, batch, BestPerSequence())
        @test size(pfm) == (4, 3)
    end

    @testset "compare and prepared profile" begin
        m1 = MatchCounter("q", Mimosa.encode_sequence("ACGT"))
        m2 = MatchCounter("t", Mimosa.encode_sequence("ACGT"))
        sequences = Mimosa.make_random_sequences(8, 60; seed=21)

        result = compare(m1, m2, sequences; metric=:co, search_range=3, window_radius=2)
        @test result isa ComparisonResult
        @test result.query == "q"
        @test result.target == "t"

        prepared = prepare_profile(m1, sequences)
        @test prepared isa PreparedProfile

        res2 = compare(prepared, m2, sequences; metric=:co, search_range=3, window_radius=2)
        @test res2.query == "q"

        targets = AbstractMotifModel[
            m2, MatchCounter("other", Mimosa.encode_sequence("ACGA"))
        ]
        results = compare(
            prepared, targets, sequences; metric=:co, search_range=3, window_radius=2
        )
        @test length(results) == 2
    end

    @testset "custom vs built-in PWM" begin
        custom = MatchCounter("custom", Mimosa.encode_sequence("ACGT"))
        examples = joinpath(dirname(dirname(@__DIR__)), "examples")
        pwm = readmodel(joinpath(examples, "pif4.meme"))
        sequences = Mimosa.make_random_sequences(6, 60; seed=33)

        r1 = compare(custom, pwm, sequences; metric=:co, search_range=3, window_radius=2)
        r2 = compare(pwm, custom, sequences; metric=:co, search_range=3, window_radius=2)
        @test r1.query == "custom"
        @test r1.target == pwm.name
        @test r2.query == pwm.name
        @test r2.target == "custom"
    end

    @testset "sites and PFM reconstruction" begin
        m = MatchCounter("downstream", Mimosa.encode_sequence("ACGT"))
        batch = Mimosa.make_random_sequences(15, 100; seed=5)

        sites = selectsites(m, batch, BestPerSequence(); strands=BothStrands())
        @test sites isa SiteCollection

        if !isempty(sites)
            pfm = reconstruct_pfm(m, batch, BestPerSequence(); pseudocount=Float32(1e-4))
            @test size(pfm) == (4, 4)
        end
    end

    @testset "compare without fingerprint" begin
        m1 = MatchCounter("q", Mimosa.encode_sequence("ACGT"))
        m2 = MatchCounter("t", Mimosa.encode_sequence("ACGT"))
        sequences = Mimosa.make_random_sequences(4, 40; seed=8)
        @test_nowarn compare(m1, m2, sequences; metric=:co, search_range=2, window_radius=1)
    end
end

end # module

@testset "Downstream contract: custom model" begin
    DownstreamCustomModel.run()
end
