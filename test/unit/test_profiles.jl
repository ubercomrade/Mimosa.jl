using Test
using Mimosa

const ExactEmpiricalLogTail = Mimosa.EmpiricalLogTail

@testset "ScoreProfile validation" begin
    @test_throws ModelFormatError ScoreProfile("bad", build_ragged([Float32[0.1, NaN32]]))
end

@testset "Hybrid normalization is the API default" begin
    profile = ScoreProfile("default", build_ragged([Float32[0.1, 0.4, 0.2]]))
    prepared = prepare_profile(profile; min_logfpr=1.0)
    @test prepared.normalization == HybridEmpiricalLogTail()
end

@testset "Hybrid normalization serial/threaded equivalence" begin
    scores = Float32[sin(index / 7) + (index % 11) / 10 for index in 1:10_000]
    strategy = HybridEmpiricalLogTail(256)
    serial = fit(strategy, scores; tail_logfpr=2.0, execution=Execution())
    threaded = fit(strategy, scores; tail_logfpr=2.0, execution=Execution(4))

    @test threaded.minimum == serial.minimum
    @test threaded.bin_width == serial.bin_width
    @test threaded.log_tail == serial.log_tail
    @test threaded.exact_tail.scores == serial.exact_tail.scores
    @test threaded.exact_tail.log_tail == serial.exact_tail.log_tail
    @test_throws ArgumentError fit(strategy, Float32[1.0, NaN32]; execution=Execution(2))
end

@testset "LogTailTable fit" begin
    # Empty input
    t = fit(ExactEmpiricalLogTail(), Float32[])
    @test length(t.scores) == 1
    @test t.scores[1] == 0.0f0
    @test t.log_tail[1] == 0.0f0

    # Single element
    t = fit(ExactEmpiricalLogTail(), Float32[5.0])
    @test t.scores == [5.0f0]
    @test t.log_tail ≈ [-0.0f0] atol = 1e-6

    # Multiple elements
    t = fit(ExactEmpiricalLogTail(), Float32[1.0, 2.0, 3.0, 1.0, 2.0, 1.0])
    # Sorted descending: [3,2,2,1,1,1]
    # Unique: [3,2,1], counts: [1,2,3]
    # Cum: [1,3,6], tail_prob: [1/6,3/6,6/6]
    # -log10: [0.778,0.301,0]
    @test t.scores == [3.0f0, 2.0f0, 1.0f0]
    @test t.log_tail[1] ≈ Float32(-log10(1.0 / 6.0)) atol = 1e-5
    @test t.log_tail[2] ≈ Float32(-log10(3.0 / 6.0)) atol = 1e-5
    @test t.log_tail[3] ≈ Float32(-log10(6.0 / 6.0)) atol = 1e-5

    input = Float32[3, 2, 2, 1, -0.0f0, 0.0f0]
    original = copy(input)
    @test fit(ExactEmpiricalLogTail(), input).scores == Float32[3, 2, 1, 0]
    @test input == original
end

@testset "LogTailTable lookup" begin
    t = fit(ExactEmpiricalLogTail(), Float32[1.0, 2.0, 3.0, 1.0, 2.0, 1.0])
    # scores descending: [3,2,1], log_tail: [0.778,0.301,0]

    # Target >= largest → index 1
    @test lookup_score(t, 3.0f0) ≈ t.log_tail[1]
    @test lookup_score(t, 5.0f0) ≈ t.log_tail[1]

    # Target <= smallest → index 3
    @test lookup_score(t, 1.0f0) ≈ t.log_tail[3]
    @test lookup_score(t, 0.0f0) ≈ t.log_tail[3]

    # Target between unique scores
    # descending scores [3,2,1]: _lower_bound_desc finds first score <= target
    @test lookup_score(t, 2.5f0) ≈ t.log_tail[2]  # 2.5 → scores[2]=2 <= 2.5
    @test lookup_score(t, 2.0f0) ≈ t.log_tail[2]  # 2.0 → scores[2]=2 <= 2.0
    @test lookup_score(t, 1.5f0) ≈ t.log_tail[3]  # 1.5 → scores[3]=1 <= 1.5
end

@testset "transform_scores" begin
    t = fit(ExactEmpiricalLogTail(), Float32[1.0, 2.0, 3.0])
    # scores: [3,2,1], log_tail: [-log10(1/3), -log10(2/3), -log10(3/3)]
    # = [0.477, 0.176, 0]

    rag = build_ragged([Float32[3.0, 2.0, 1.0], Float32[2.0, 1.0]])
    transformed = transform_scores(t, rag)

    @test nrows(transformed) == 2
    @test rowlength(transformed, 1) == 3
    @test rowlength(transformed, 2) == 2
    @test row(transformed, 1) ≈ [t.log_tail[1], t.log_tail[2], t.log_tail[3]]
    @test row(transformed, 2) ≈ [t.log_tail[2], t.log_tail[3]]
end

@testset "sorted transform matches scalar lookup" begin
    table = fit(ExactEmpiricalLogTail(), Float32[4, 4, 3, 1, -2, -2])
    rag = build_ragged([
        Float32[8, 4, 3.5, 3, 2, 1, 0, -2, -3], Float32[4, 4, -2], Float32[]
    ])
    expected = Float32[lookup_score(table, score) for score in rag.data]

    transformed = transform_scores(table, rag)
    @test transformed.data == expected
    @test transformed.offsets == rag.offsets

    if Threads.nthreads() > 1
        threaded = transform_scores(table, rag; execution=Execution(2))
        @test threaded == transformed
    end
end

@testset "normalize_bundle" begin
    # Build a simple strand bundle
    rag = build_ragged([Float32[3.0, 1.0], Float32[2.0, 2.0]])
    bundle = StrandPair(rag, rag)

    flat = flatten_bundle(bundle)
    @test length(flat) == 8  # 4 elements × 2 strands

    t = fit(ExactEmpiricalLogTail(), flat)
    normed = normalize_bundle(t, bundle)

    # Both strands should be transformed
    @test nrows(normed.forward) == 2
    @test nrows(normed.reverse) == 2
    @test row(normed.forward, 1) ≈ row(normed.reverse, 1)
end

@testset "fused two-strand empirical normalization" begin
    forward = build_ragged([Float32[3, 1, 3], Float32[], Float32[2, 0]])
    reverse = build_ragged([Float32[0, 2, 3], Float32[1], Float32[]])
    bundle = StrandPair(forward, reverse)

    fused_table, fused = Mimosa._fit_normalize_empirical(bundle)
    reference_table = fit(ExactEmpiricalLogTail(), flatten_bundle(bundle))
    reference = normalize_bundle(reference_table, bundle)

    @test fused_table.scores == reference_table.scores
    @test fused_table.log_tail == reference_table.log_tail
    @test fused.forward.offsets == forward.offsets
    @test fused.reverse.offsets == reverse.offsets
    @test fused.forward.data == reference.forward.data
    @test fused.reverse.data == reference.reverse.data

    _, symmetric = Mimosa._fit_normalize_empirical(StrandPair(forward, forward))
    @test symmetric.forward === symmetric.reverse
end

@testset "best-anchor scalar alignment matches CSR alignment" begin
    query = build_ragged([Float32[0.1, 0.8, 0.2, 0.4, 0.3], Float32[]])
    target = build_ragged([Float32[0.2, 0.7, 0.1, 0.5, 0.3], Float32[]])
    query_anchors = Mimosa._collect_both_anchors(StrandPair(query, query), 0.0f0)[1]
    target_anchors = Mimosa._collect_both_anchors(StrandPair(target, target), 0.0f0)[1]
    metrics = (
        OverlapCoefficient(),
        OverlapCoefficientRowwise(),
        DiceSimilarity(),
        DiceSimilarityRowwise(),
        CosineSimilarityProfile(),
    )

    for metric in metrics
        for shift in -2:2
            scratch = Mimosa.CandidateScratch(5)
            reference = Mimosa._score_shift!(
                scratch, query, target, query_anchors, target_anchors, shift, 1, 1, metric
            )
            optimized = Mimosa._score_shift_best!(
                query, target, query_anchors, target_anchors, shift, 1, 1, metric
            )
            @test optimized == reference
        end
    end
end

@testset "AnchorCSR" begin
    # Build anchors: row 1 has positions [3, 5], row 2 has position [1]
    rows = [1, 1, 2]
    positions = [3, 5, 1]
    csr = build_anchor_csr(rows, positions, 3)

    @test csr.offsets == [1, 3, 4, 4]  # row 1: [1:2], row 2: [3:3], row 3: empty
    @test csr.positions == [3, 5, 1]  # stable sort preserves order

    # Empty case
    csr_empty = build_anchor_csr(Int[], Int[], 2)
    @test isempty(csr_empty)
    @test csr_empty.offsets == [1, 1, 1]

    @test_throws ArgumentError AnchorCSR([1], [2, 2])
    @test_throws ArgumentError AnchorCSR([1], [1, 3])
    @test_throws ArgumentError AnchorCSR([1], [1, 0])
    @test_throws ArgumentError build_anchor_csr([1], Int[], 1)
    @test_throws ArgumentError build_anchor_csr([0], [1], 1)
    @test_throws ArgumentError build_anchor_csr([2], [1], 1)
    @test_throws ArgumentError build_anchor_csr([1], [0], 1)
    @test_throws ArgumentError build_anchor_csr(Int[], Int[], -1)
end

@testset "collect_best_anchors" begin
    rag = build_ragged([
        Float32[0.1, 0.5, 0.3],  # best at position 2
        Float32[0.7, 0.2, 0.9, 0.1],  # best at position 3
        Float32[],  # empty row
    ])

    rows, positions = collect_best_anchors(rag)
    @test rows == [1, 2]
    @test positions == [2, 3]
end

@testset "collect_threshold_anchors" begin
    rag = build_ragged([Float32[0.1, 0.5, 0.3], Float32[0.7, 0.2, 0.9, 0.1]])

    rows, positions = collect_threshold_anchors(rag, 0.4f0)
    @test rows == [1, 2, 2]
    @test positions == [2, 1, 3]
end

@testset "anchor CSR serial/threaded equivalence" begin
    scores = build_ragged([Float32[0.1, 0.9, 0.9], Float32[], Float32[0.8, 0.2, 0.7, 0.8]])
    for threshold in (0.0f0, 0.75f0)
        serial = Mimosa.collect_anchor_csr(scores, threshold)
        threaded = Mimosa.collect_anchor_csr(scores, threshold; execution=Execution(4))
        @test threaded.positions == serial.positions
        @test threaded.offsets == serial.offsets
    end
end

@testset "ProfileConfig defaults" begin
    config = ProfileConfig()
    @test config.metric isa OverlapCoefficient
    @test config.search_range == 10
    @test config.window_radius == 10
    @test config.realign_window == 3
    @test config.min_logfpr == 0.0f0

    config2 = ProfileConfig(metric=DiceSimilarity(), search_range=5, window_radius=3)
    @test config2.metric isa DiceSimilarity
    @test config2.search_range == 5
    @test config2.window_radius == 3
end

@testset "Profile orientation priority is explicit" begin
    @test [p[1] for p in Mimosa.PROFILE_ORIENTATION_PAIRS] == ["++", "+-", "-+", "--"]
    @test Mimosa.PROFILE_ORIENTATION_RANK["+-"] < Mimosa.PROFILE_ORIENTATION_RANK["--"]
    @test Mimosa.PROFILE_ORIENTATION_RANK["-+"] < Mimosa.PROFILE_ORIENTATION_RANK["--"]
    @test_throws ArgumentError ProfileConfig(search_range=-1)
    @test_throws ArgumentError ProfileConfig(min_logfpr=NaN)
end

@testset "PreparedProfile and one-to-many" begin
    # Create two ScoreProfiles
    sp1 = ScoreProfile(
        "query", build_ragged([Float32[0.1, 0.5, 0.3, 0.8, 0.2, 0.6, 0.1, 0.9, 0.3, 0.4]])
    )
    sp2 = ScoreProfile(
        "target", build_ragged([Float32[0.2, 0.4, 0.3, 0.7, 0.3, 0.5, 0.2, 0.8, 0.4, 0.3]])
    )
    sp3 = ScoreProfile(
        "target2", build_ragged([Float32[0.3, 0.1, 0.9, 0.2, 0.8, 0.1, 0.5, 0.3, 0.7, 0.2]])
    )

    @test sp1 isa AbstractProfileSource
    @test !(sp1 isa AbstractMotifModel)
    @test !applicable(motif_length, sp1)

    # Prepare the query
    prepared = prepare_profile(sp1)
    @test prepared.name == "query"
    @test prepared.bundle isa StrandPair
    @test prepared.anchors isa Tuple{AnchorCSR,AnchorCSR}

    # Compare prepared vs ScoreProfile should match direct compare
    direct = compare(sp1, sp2; metric=:co, search_range=3, window_radius=2)
    prepared_result = compare(prepared, sp2; metric=:co, search_range=3, window_radius=2)
    @test prepared_result.score ≈ direct.score atol = 1e-5
    @test prepared_result.offset == direct.offset
    @test prepared_result.orientation == direct.orientation

    # One-to-many comparison
    results = compare(prepared, [sp2, sp3]; metric=:co, search_range=3, window_radius=2)
    @test length(results) == 2
    @test results[1].score ≈ prepared_result.score atol = 1e-5
    @test results[1].offset == prepared_result.offset
    @test results[2].query == "query"
    @test results[2].target == "target2"

    thresholded = prepare_profile(sp1; min_logfpr=0.25)
    threshold_serial = compare(thresholded, sp2; search_range=3, window_radius=2)
    threshold_threaded = compare(
        thresholded, sp2; search_range=3, window_radius=2, execution=Execution(4)
    )
    @test threshold_serial.score isa Float32
    @test threshold_threaded == threshold_serial
    @test_throws ArgumentError compare(
        thresholded, sp2; search_range=3, window_radius=2, min_logfpr=0.0
    )
    @test_throws ArgumentError compare(
        thresholded, [sp2]; search_range=3, window_radius=2, min_logfpr=0.0
    )

    # Determinism: repeated calls give same results
    results2 = compare(prepared, [sp2, sp3]; metric=:co, search_range=3, window_radius=2)
    @test results[1].score ≈ results2[1].score atol = 1e-6
    @test results[2].score ≈ results2[2].score atol = 1e-6
end

@testset "motif model one-to-many profile comparison" begin
    weights = Float32[
        0.5 -0.2 0.1
        -0.1 0.6 -0.3
        0.2 -0.1 0.7
        -0.4 0.1 -0.2
        -0.5 -0.5 -0.5
    ]
    background = (0.25f0, 0.25f0, 0.25f0, 0.25f0)
    query = PWM("query", weights, background)
    targets = [
        PWM("target1", weights .+ 0.05f0, background),
        PWM("target2", reverse(weights; dims=2), background),
    ]
    batch = EncodedSequenceBatch(
        UInt8[0, 1, 2, 3, 0, 1, 2, 3, 3, 2, 1, 0, 3, 2, 1, 0], [1, 9, 17]
    )

    serial = compare(
        query,
        targets,
        batch;
        metric=:co,
        search_range=2,
        window_radius=1,
        realign_window=1,
        execution=Execution(),
    )
    prepared = compare(
        prepare_profile(query, batch),
        targets,
        batch;
        metric=:co,
        search_range=2,
        window_radius=1,
        realign_window=1,
        execution=Execution(),
    )
    @test [r.target for r in serial] == ["target1", "target2"]
    @test prepared == serial

    threaded = compare(
        query,
        targets,
        batch;
        metric=:co,
        search_range=2,
        window_radius=1,
        realign_window=1,
        execution=Execution(4),
    )
    @test threaded == serial

    scan_threaded = compare(
        query,
        targets[1],
        batch;
        metric=:co,
        search_range=2,
        window_radius=1,
        realign_window=1,
        execution=Execution(2),
    )
    @test scan_threaded == serial[1]
end
