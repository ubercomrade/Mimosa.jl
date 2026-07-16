using Test
using Mimosa

const _NULL_BG = (0.25f0, 0.25f0, 0.25f0, 0.25f0)

@testset "profile-only null configuration" begin
    config = NullBuildConfig(metric=:co, min_null_targets=1)
    @test config.metric isa OverlapCoefficient
    @test config.min_null_targets == 1
    @test_throws ArgumentError NullBuildConfig(metric=:not_a_profile_metric)
    @test_throws ArgumentError NullBuildConfig(metric=:co, min_null_targets=0)
end

@testset "profile null build accepts precomputed profile sources" begin
    models = [
        ScoreProfile("p1", build_ragged([Float32[0.1, 0.8, 0.2, 0.7]])),
        ScoreProfile("p2", build_ragged([Float32[0.2, 0.7, 0.3, 0.6]])),
    ]
    relations = GroupRelations(
        Dict("p1" => "A", "p2" => "B"), Dict("p1" => Set(["p2"]), "p2" => Set(["p1"]))
    )
    sequences = EncodedSequenceBatch([encode_sequence("ACGT")])

    result = build_null(models, relations; sequences=sequences, metric=:co)
    @test result.total_comparisons == 2
    @test result.distribution.model_collection_fingerprint ==
        model_collection_fingerprint(models)
end

@testset "null annotation validates metric" begin
    skipped = NamedTuple{(:query, :reason),Tuple{String,String}}[]
    fit = GEVFit(0.0, 0.0, 1.0, true, 1, -1.0)
    dist = NullDistribution(
        "profile",
        "co",
        fit,
        Float64[0, 1, 2],
        NullPair[],
        3,
        1,
        skipped,
        nothing,
        nothing,
        "none",
        "none",
    )
    wrong_metric = ComparisonResult("q", "t", 0.5f0, 0, "++", "dice", 1)
    @test_throws ArgumentError annotate_results([wrong_metric], dist)
end

@testset "profile null build" begin
    w1 = Float32[0.8 0.1; 0.1 0.8; 0.05 0.05; 0.05 0.05; 0.0 0.0]
    w2 = Float32[0.1 0.8; 0.8 0.1; 0.05 0.05; 0.05 0.05; 0.0 0.0]
    models = [PWM("m1", w1, _NULL_BG), PWM("m2", w2, _NULL_BG)]
    relations = GroupRelations(
        Dict("m1" => "A", "m2" => "B"), Dict("m1" => Set(["m2"]), "m2" => Set(["m1"]))
    )
    sequences = EncodedSequenceBatch([
        encode_sequence("ACGTACGT"), encode_sequence("TGCATGCA")
    ])
    result = build_null(models, relations; sequences=sequences, metric=:co)
    @test result.distribution.strategy == "profile"
    @test result.distribution.metric == "co"
    @test result.distribution.sequence_fingerprint == sequence_fingerprint(sequences)
    @test result.total_comparisons == 2
    duplicate = [models[1], PWM("m1", w2, _NULL_BG)]
    @test_throws ArgumentError build_null(duplicate, relations; sequences=sequences)
end
