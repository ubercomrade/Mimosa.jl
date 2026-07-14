using Mimosa
using Printf

function median_elapsed(times::Vector{Float64})
    ordered = sort(times)
    midpoint = length(ordered) ÷ 2
    return if isodd(length(ordered))
        ordered[midpoint + 1]
    else
        (ordered[midpoint] + ordered[midpoint + 1]) / 2
    end
end

const PROFILE_KWARGS = (; metric=:co, search_range=10, window_radius=5, realign_window=3)

function timed_workload(query, targets, sequences, execution)
    return compare(query, targets, sequences; execution, PROFILE_KWARGS...)
end

function stage_stats(stage, f, reps)
    f()
    times = Vector{Float64}(undef, reps)
    for i in eachindex(times)
        GC.gc()
        times[i] = @elapsed f()
    end
    ordered = sort(times)
    midpoint = length(ordered) ÷ 2
    median = if isodd(length(ordered))
        ordered[midpoint + 1]
    else
        (ordered[midpoint] + ordered[midpoint + 1]) / 2
    end
    @printf(
        "STAGE language=julia name=%s median_s=%.9f min_s=%.9f\n",
        stage,
        median,
        first(ordered)
    )
end

function main(args)
    length(args) == 4 || error("usage: cross_language_profile.jl FASTA MEME THREADS REPS")
    fasta, meme = args[1], args[2]
    threads = parse(Int, args[3])
    reps = parse(Int, args[4])
    threads <= Threads.nthreads() ||
        error("requested $threads threads, runtime has $(Threads.nthreads())")
    execution = threads == 1 ? SerialExecution() : ThreadedExecution(threads)

    sequences, _ = readsequences(fasta)
    models = [readmodel(meme; index=i) for i in 0:50]
    query, targets = first(models), models[2:end]

    raw_query = scan(query, sequences; strands=BothStrands(), execution=SerialExecution())
    raw_targets = [
        scan(target, sequences; strands=BothStrands(), execution=SerialExecution()) for
        target in targets
    ]
    _, normalized_query = Mimosa._fit_transform_empirical(raw_query)
    normalized_targets = [Mimosa._fit_transform_empirical(raw)[2] for raw in raw_targets]
    prepared_query = prepare_profile(query, sequences; execution=SerialExecution())
    prepared_targets = [
        PreparedProfile(
            target.name,
            normalized_targets[i],
            Mimosa._collect_both_anchors(normalized_targets[i], 0.0f0),
            0.0f0,
        ) for (i, target) in enumerate(targets)
    ]
    stage_stats(
        "query_scan", () -> scan(query, sequences; strands=BothStrands(), execution), reps
    )
    stage_stats(
        "query_normalization", () -> Mimosa._fit_transform_empirical(raw_query), reps
    )
    stage_stats(
        "target_scan",
        () -> [
            scan(target, sequences; strands=BothStrands(), execution=SerialExecution())
            for target in targets
        ],
        reps,
    )
    stage_stats(
        "target_normalization",
        () -> [Mimosa._fit_transform_empirical(raw)[2] for raw in raw_targets],
        reps,
    )
    stage_stats(
        "anchor_collection",
        () -> [
            Mimosa._collect_both_anchors(normalized, 0.0f0) for
            normalized in [normalized_query; normalized_targets...]
        ],
        reps,
    )
    config = ProfileConfig(;
        metric=OverlapCoefficient(),
        search_range=10,
        window_radius=5,
        realign_window=3,
        min_logfpr=0.0f0,
    )
    stage_stats(
        "alignment_1v1",
        () -> profile_compare(
            prepared_query.bundle,
            prepared_query.anchors,
            prepared_targets[1].bundle,
            prepared_targets[1].anchors,
            config,
        ),
        reps,
    )
    stage_stats(
        "prepared_1v50",
        () -> [
            compare(prepared_query, target; PROFILE_KWARGS...) for
            target in prepared_targets
        ],
        reps,
    )

    length(timed_workload(query, targets, sequences, execution)) == 50 ||
        error("warm-up failed")
    times = Vector{Float64}(undef, reps)
    for i in eachindex(times)
        GC.gc()
        times[i] = @elapsed timed_workload(query, targets, sequences, execution)
    end
    @printf(
        "RESULT language=julia threads=%d median_s=%.9f min_s=%.9f julia_version=%s\n",
        threads,
        median_elapsed(times),
        minimum(times),
        VERSION,
    )
end

main(ARGS)
