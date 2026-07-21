#!/usr/bin/env julia
# Benchmark Hybrid serial/threaded scaling and the internal Exact reference:
#
#   1. Hybrid serial and inner-threaded
#      - HybridEmpiricalLogTail normalization;
#      - targets are processed serially;
#      - computational kernels inside each target use ThreadedExecution.
#
#   2. Exact inner-threaded
#      - EmpiricalLogTail normalization;
#      - targets are processed serially;
#      - computational kernels inside each target use ThreadedExecution.
#
# Environment overrides:
#   MIMOSA_BENCH_N_SEQUENCES, MIMOSA_BENCH_SEQ_LENGTH,
#   MIMOSA_BENCH_N_TARGETS, MIMOSA_BENCH_PWM_WIDTH, MIMOSA_BENCH_N_REPS,
#   MIMOSA_BENCH_HYBRID_BINS, MIMOSA_BENCH_MIN_LOGFPR, MIMOSA_BENCH_SEED.

using Dates
using Mimosa
using Printf
using Random

# ── Helpers ──────────────────────────────────────────────────────────────────

function make_pwm(width::Int, seed::Int=42)
    rng = Random.MersenneTwister(seed)
    weights = Matrix{Float32}(undef, 5, width)
    for col in 1:width
        for row in 1:4
            weights[row, col] = Float32(randn(rng) * 0.5)
        end
        weights[5, col] = minimum(@view(weights[1:4, col])) - 0.1f0
    end
    background = (0.25f0, 0.25f0, 0.25f0, 0.25f0)
    return PWM("pwm_w$(width)_s$seed", weights, background)
end

function elapsed(f::F) where {F}
    start = Base.time_ns()
    value = f()
    return (nanoseconds=Base.time_ns() - start, value=value)
end

function bench_repeat(f::F, repetitions::Int) where {F}
    f() # compile and warm up outside measured samples
    times = Vector{Int}(undef, repetitions)
    for index in eachindex(times)
        GC.gc()
        times[index] = elapsed(f).nanoseconds
    end
    return times
end

function stats(times::Vector{Int})
    sorted = sort(times)
    median = sorted[cld(length(sorted), 2)]
    return (
        min_ns=first(sorted),
        median_ns=median,
        max_ns=last(sorted),
        mean_ns=sum(times) / length(times),
    )
end

fmt_ms(nanoseconds::Real) = @sprintf("%.3f", nanoseconds / 1.0e6)

struct BenchmarkCase{N<:AbstractNormalizationStrategy,E<:ExecutionPolicy}
    label::String
    normalization::N
    execution::E
end

function show_case(case::BenchmarkCase)
    println("  $(case.label)")
    println("    normalization : $(normalization_fingerprint(case.normalization))")
    println("    targets       : serial")
    return println("    kernels       : $(case.execution)")
end

function benchmark_case(
    case::BenchmarkCase,
    query_model::AbstractMotifModel,
    target_models::AbstractVector{<:AbstractMotifModel},
    sequences::EncodedSequenceBatch,
    repetitions::Int,
    ;
    background::Union{EncodedSequenceBatch,Nothing}=nothing,
)
    println("\n[$(case.label)] Preparing query...")
    prepare_query() = prepare_profile(
        query_model,
        sequences;
        background=background,
        min_logfpr=MIN_LOGFPR,
        normalization=case.normalization,
        execution=case.execution,
    )
    query_times = bench_repeat(prepare_query, repetitions)
    query_stats = stats(query_times)
    prepared_query = prepare_query()
    println(
        "  query: min=$(fmt_ms(query_stats.min_ns)) ms  " *
        "median=$(fmt_ms(query_stats.median_ns)) ms  " *
        "max=$(fmt_ms(query_stats.max_ns)) ms",
    )

    compare_one() = compare(
        prepared_query,
        target_models[1],
        sequences;
        background=background,
        execution=case.execution,
        metric=:co,
        search_range=10,
        window_radius=5,
    )
    one_times = bench_repeat(compare_one, repetitions)
    one_stats = stats(one_times)
    println(
        "  1-vs-1: min=$(fmt_ms(one_stats.min_ns)) ms  " *
        "median=$(fmt_ms(one_stats.median_ns)) ms  " *
        "max=$(fmt_ms(one_stats.max_ns)) ms",
    )

    compare_many() = compare(
        prepared_query,
        target_models,
        sequences;
        background=background,
        execution=case.execution,
        metric=:co,
        search_range=10,
        window_radius=5,
    )
    many_times = bench_repeat(compare_many, repetitions)
    many_stats = stats(many_times)
    results = compare_many()
    per_target_ms = many_stats.median_ns / length(target_models) / 1.0e6
    throughput = length(target_models) / (many_stats.median_ns / 1.0e9)
    println(
        "  1-vs-$(length(target_models)): min=$(fmt_ms(many_stats.min_ns)) ms  " *
        "median=$(fmt_ms(many_stats.median_ns)) ms  " *
        "max=$(fmt_ms(many_stats.max_ns)) ms",
    )
    println(
        @sprintf(
            "  per target: %.3f ms; throughput: %.1f targets/s", per_target_ms, throughput
        )
    )

    return (case=case, query=query_stats, one=one_stats, many=many_stats, results=results)
end

function compare_outputs(exact_results, hybrid_results)
    length(exact_results) == length(hybrid_results) ||
        error("benchmark strategies returned different result counts")
    score_deltas = abs.(
        getproperty.(exact_results, :score) .- getproperty.(hybrid_results, :score)
    )
    offset_changes = count(
        index -> exact_results[index].offset != hybrid_results[index].offset,
        eachindex(exact_results),
    )
    orientation_changes = count(
        index -> exact_results[index].orientation != hybrid_results[index].orientation,
        eachindex(exact_results),
    )
    site_changes = count(
        index -> exact_results[index].n_sites != hybrid_results[index].n_sites,
        eachindex(exact_results),
    )
    return (
        mean_score_delta=sum(score_deltas) / length(score_deltas),
        max_score_delta=maximum(score_deltas),
        offset_changes=offset_changes,
        orientation_changes=orientation_changes,
        site_changes=site_changes,
    )
end

# ── Configuration ────────────────────────────────────────────────────────────

const SEQ_LENGTH = parse(Int, get(ENV, "MIMOSA_BENCH_SEQ_LENGTH", "100"))
const N_SEQUENCES = parse(Int, get(ENV, "MIMOSA_BENCH_N_SEQUENCES", "10000"))
const N_TARGETS = parse(Int, get(ENV, "MIMOSA_BENCH_N_TARGETS", "50"))
const PWM_WIDTH = parse(Int, get(ENV, "MIMOSA_BENCH_PWM_WIDTH", "15"))
const N_REPS = parse(Int, get(ENV, "MIMOSA_BENCH_N_REPS", "5"))
const HYBRID_BINS = parse(Int, get(ENV, "MIMOSA_BENCH_HYBRID_BINS", "65536"))
const MIN_LOGFPR = parse(Float32, get(ENV, "MIMOSA_BENCH_MIN_LOGFPR", "3.0"))
const SEED = parse(Int, get(ENV, "MIMOSA_BENCH_SEED", "12345"))

# ── Main benchmark ───────────────────────────────────────────────────────────

function main()
    N_SEQUENCES > 0 || error("MIMOSA_BENCH_N_SEQUENCES must be positive")
    N_TARGETS > 0 || error("MIMOSA_BENCH_N_TARGETS must be positive")
    N_REPS > 0 || error("MIMOSA_BENCH_N_REPS must be positive")

    threaded = ThreadedExecution()
    hybrid_serial_case = BenchmarkCase(
        "Hybrid serial", HybridEmpiricalLogTail(HYBRID_BINS), SerialExecution()
    )
    hybrid_threaded_case = BenchmarkCase(
        "Hybrid inner-threaded", HybridEmpiricalLogTail(HYBRID_BINS), threaded
    )
    exact_case = BenchmarkCase("Exact inner-threaded", Mimosa.EmpiricalLogTail(), threaded)

    println("=" ^ 76)
    println("  Mimosa.jl — 1-vs-$N_TARGETS normalization/execution benchmark")
    println("=" ^ 76)
    println("  Sequence length    : $SEQ_LENGTH")
    println("  Number of sequences: $N_SEQUENCES")
    println("  Number of targets  : $N_TARGETS")
    println("  PWM width          : $PWM_WIDTH")
    println("  Hybrid bins        : $HYBRID_BINS")
    println("  Minimum log FPR    : $MIN_LOGFPR")
    println("  Repetitions        : $N_REPS")
    println("  Seed               : $SEED")
    println("  Julia threads      : $(Threads.nthreads())")
    println("  Date               : $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
    println("\nExecution cases:")
    show_case(hybrid_serial_case)
    show_case(hybrid_threaded_case)
    show_case(exact_case)
    println("=" ^ 76)

    generated = elapsed() do
        return make_random_sequences(N_SEQUENCES, SEQ_LENGTH; seed=SEED)
    end
    sequences = generated.value
    println(
        "\nGenerated $N_SEQUENCES × $SEQ_LENGTH bp in " *
        "$(fmt_ms(generated.nanoseconds)) ms",
    )

    query_model = make_pwm(PWM_WIDTH, SEED)
    target_models = [make_pwm(PWM_WIDTH, SEED + 100 + index) for index in 1:N_TARGETS]

    # Every case performs its own warm-up before measured samples.
    hybrid_serial = benchmark_case(
        hybrid_serial_case, query_model, target_models, sequences, N_REPS
    )
    hybrid_threaded = benchmark_case(
        hybrid_threaded_case, query_model, target_models, sequences, N_REPS
    )
    exact = benchmark_case(exact_case, query_model, target_models, sequences, N_REPS)
    hybrid_serial.results == hybrid_threaded.results ||
        error("serial and threaded Hybrid results differ")
    differences = compare_outputs(exact.results, hybrid_threaded.results)

    hybrid_serial_ns = Float64(hybrid_serial.many.median_ns)
    hybrid_threaded_ns = Float64(hybrid_threaded.many.median_ns)
    exact_ns = Float64(exact.many.median_ns)
    println("\n" * "=" ^ 76)
    println("  SUMMARY")
    println("=" ^ 76)
    println("  Hybrid serial          1-vs-$N_TARGETS: $(fmt_ms(hybrid_serial_ns)) ms")
    println(
        "  Hybrid inner-threaded  1-vs-$N_TARGETS: $(fmt_ms(hybrid_threaded_ns)) ms " *
        "($(normalization_fingerprint(hybrid_threaded_case.normalization)))",
    )
    println(
        "  Exact inner-threaded   1-vs-$N_TARGETS: $(fmt_ms(exact_ns)) ms " *
        "($(normalization_fingerprint(exact_case.normalization)))",
    )
    println(
        @sprintf("  Hybrid kernel speedup: %.3fx", hybrid_serial_ns / hybrid_threaded_ns,)
    )
    println(@sprintf("  Hybrid / Exact time ratio: %.3fx", hybrid_threaded_ns / exact_ns))
    println("\n  Result differences (Hybrid versus Exact):")
    println(@sprintf("    mean |Δscore| : %.7f", differences.mean_score_delta))
    println(@sprintf("    max  |Δscore| : %.7f", differences.max_score_delta))
    println("    changed offsets     : $(differences.offset_changes)/$N_TARGETS")
    println("    changed orientations: $(differences.orientation_changes)/$N_TARGETS")
    println("    changed site counts : $(differences.site_changes)/$N_TARGETS")

    println("\n  Sample results (first $(min(5, N_TARGETS)) targets):")
    for index in 1:min(5, N_TARGETS)
        exact_result = exact.results[index]
        hybrid_result = hybrid_threaded.results[index]
        println(
            "    $(exact_result.target): " * @sprintf(
                "exact=%.4f hybrid=%.4f Δ=%+.5f",
                exact_result.score,
                hybrid_result.score,
                hybrid_result.score - exact_result.score,
            ),
        )
    end
    println("=" ^ 76)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
