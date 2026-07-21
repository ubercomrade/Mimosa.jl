using Test
using Mimosa

const _NULL_BG = (0.25f0, 0.25f0, 0.25f0, 0.25f0)

function _null_test_models()
    w1 = Float32[
        0.8 0.1 0.2
        0.1 0.8 0.3
        0.05 0.05 0.4
        0.05 0.05 0.1
        0.05 0.05 0.1
    ]
    w2 = Float32[
        0.1 0.8 0.4
        0.8 0.1 0.2
        0.05 0.05 0.3
        0.05 0.05 0.1
        0.05 0.05 0.1
    ]
    return [PWM("m1", w1, _NULL_BG), PWM("m2", w2, _NULL_BG)]
end

@testset "random-pair null configuration" begin
    config = NullBuildConfig(metric=:co, n_samples=17, shuffle=true, seed=42)
    @test config.metric isa OverlapCoefficient
    @test config.n_samples == 17
    @test config.shuffle
    @test config.seed == 42
    @test_throws ArgumentError NullBuildConfig(metric=:not_a_profile_metric)
    @test_throws ArgumentError NullBuildConfig(n_samples=0)
    @test_throws ArgumentError NullBuildConfig(seed=-1)
end

@testset "PWM null shuffling" begin
    model = first(_null_test_models())
    shuffled = Mimosa._shuffle_null_model(model, UInt64(123))

    @test shuffled !== model
    @test modelname(shuffled) == modelname(model)
    @test shuffled.background == model.background
    @test size(shuffled) == size(model)
    @test shuffled.representation != model.representation
    @test all(
        shuffled.representation[5, column] ==
        minimum(@view shuffled.representation[1:4, column]) for
        column in axes(shuffled.representation, 2)
    )

    original_columns = sort([sort(model.representation[1:4, column]) for column in 1:3])
    shuffled_columns = sort([sort(shuffled.representation[1:4, column]) for column in 1:3])
    @test shuffled_columns == original_columns
    @test Mimosa._shuffle_null_model(model, UInt64(123)) == shuffled

    bamm = BaMM("bamm", zeros(Float32, 5, 3), 0, 3)
    @test Mimosa._shuffle_null_model(bamm, UInt64(123)) === bamm
end

@testset "random-pair profile null build" begin
    models = _null_test_models()
    sequences = EncodedSequenceBatch([
        encode_sequence("ACGTACGT"), encode_sequence("TGCATGCA")
    ])
    progress_events = NamedTuple[]
    result = build_null(
        models;
        sequences=sequences,
        metric=:co,
        n_samples=12,
        shuffle=true,
        seed=9,
        on_progress=event -> push!(progress_events, event),
    )

    dist = result.distribution
    @test result.total_comparisons == 12
    @test dist.n_null == 12
    @test dist.n_models == 2
    @test dist.model_type == "pwm"
    @test dist.shuffle
    @test dist.seed == 9
    @test dist.sampling_version == "random-ordered-pairs-v1"
    @test length(dist.pairs) == 12
    @test all(pair.query != pair.target for pair in dist.pairs)
    @test dist.sequence_fingerprint == sequence_fingerprint(sequences)
    @test dist.model_collection_fingerprint == model_collection_fingerprint(models)
    @test all(event.stage === :null for event in progress_events)
    @test [event.current for event in progress_events] == collect(0:12)
    @test all(event.total == 12 for event in progress_events)

    staged_events = NamedTuple[]
    build_null(
        models;
        sequences=sequences,
        metric=:co,
        n_samples=2,
        shuffle=false,
        seed=9,
        on_progress=event -> push!(staged_events, event),
    )
    @test [(event.stage, event.current, event.total) for event in staged_events] == [
        (:prepare, 0, 2),
        (:prepare, 1, 2),
        (:prepare, 2, 2),
        (:null, 0, 2),
        (:null, 1, 2),
        (:null, 2, 2),
    ]

    repeated = build_null(
        models; sequences=sequences, metric=:co, n_samples=12, shuffle=true, seed=9
    )
    @test repeated.distribution.raw_scores == dist.raw_scores
    @test repeated.distribution.pairs == dist.pairs

    scan_threaded = build_null(
        models;
        sequences=sequences,
        metric=:co,
        n_samples=12,
        shuffle=true,
        seed=9,
        execution=Execution(2),
    )
    @test scan_threaded.distribution.raw_scores == dist.raw_scores
    @test scan_threaded.distribution.pairs == dist.pairs
    @test_throws ArgumentError build_null(models[1:1]; sequences=sequences, n_samples=2)
    duplicate = [models[1], PWM("m1", models[2].representation, _NULL_BG)]
    @test_throws ArgumentError build_null(duplicate; sequences=sequences, n_samples=2)
end

@testset "null annotation validates metric" begin
    fit = GEVFit(0.0, 0.0, 1.0, true, 1, -1.0)
    scores = Float64[0, 1, 2]
    dist = NullDistribution(
        "profile",
        "co",
        fit,
        scores,
        [NullPair("q", "t", score) for score in scores],
        3,
        2,
        "pwm",
        false,
        127,
        "random-ordered-pairs-v1",
        nothing,
        "none",
        "none",
    )
    wrong_metric = ComparisonResult("q", "t", 0.5f0, 0, "++", "dice", 1)
    @test_throws ArgumentError annotate_results([wrong_metric], dist)
end
