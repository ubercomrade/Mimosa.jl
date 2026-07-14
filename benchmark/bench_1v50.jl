#!/usr/bin/env julia
# Benchmark: 1-vs-50 profile comparison on random sequences.
#
# Parameters (per user request):
#   sequence length  = 100
#   number of sequences = 10000
#   comparison = 1 query vs 50 targets
#
# Measures:
#   1. Sequence generation time
#   2. Query model scan time (against 10000 x 100 sequences)
#   3. Query profile preparation time
#   4. Single target scan time
#   5. 1-vs-1 comparison against a raw ScoreProfile
#   6. 1-vs-50 comparison against raw ScoreProfiles
#   7. Full end-to-end 1-vs-50 (scan each target + compare, query pre-prepared)
#   8. Throughput (comparisons/second)

using Random
using Printf
using Dates
using Mimosa

# ── Helpers ──────────────────────────────────────────────────────────────────

function make_pwm(width::Int, seed::Int=42)
    rng = Random.MersenneTwister(seed)
    weights = Matrix{Float32}(undef, 5, width)
    for col in 1:width
        for row in 1:4
            weights[row, col] = Float32(randn(rng) * 0.5)
        end
        weights[5, col] = minimum(@view weights[1:4, col]) - Float32(0.1)
    end
    bg = (Float32(0.25), Float32(0.25), Float32(0.25), Float32(0.25))
    return PWM("pwm_w$(width)_s$seed", weights, bg)
end

function elapsed(f::Function)
    t0 = Base.time_ns()
    f()
    return Base.time_ns() - t0
end

function bench_repeat(f::Function, n_reps::Int)
    times = Vector{Int}(undef, n_reps)
    for i in 1:n_reps
        GC.gc()
        times[i] = elapsed(f)
    end
    return times
end

function stats(times::Vector{Int})
    sorted = sort(times)
    med = sorted[div(length(sorted) + 1, 2)]
    mn = sorted[1]
    mx = sorted[end]
    mean = sum(times) / length(times)
    return (min_ns=mn, median_ns=med, max_ns=mx, mean_ns=mean)
end

function fmt_ms(ns::Float64)
    return @sprintf("%.3f", ns / 1e6)
end

# ── Configuration ─────────────────────────────────────────────────────────────

const SEQ_LENGTH = 100
const N_SEQUENCES = 10000
const N_TARGETS = 50
const PWM_WIDTH = 15
const N_REPS = 5
const SEED = 12345

# ── Main benchmark ─────────────────────────────────────────────────────────────

function main()
    execution = Threads.nthreads() == 1 ? SerialExecution() : ThreadedExecution()
    println("=" ^ 72)
    println("  Mimosa.jl — 1-vs-50 Profile Comparison Benchmark")
    println("=" ^ 72)
    println("  Sequence length    : $SEQ_LENGTH")
    println("  Number of sequences: $N_SEQUENCES")
    println("  Number of targets  : $N_TARGETS")
    println("  PWM width          : $PWM_WIDTH")
    println("  Repetitions       : $N_REPS")
    println("  Seed               : $SEED")
    println("  Julia threads      : $(Threads.nthreads())")
    println("  Execution policy   : $execution")
    println("  Date               : $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
    println("=" ^ 72)

    # ── 1. Generate random sequences ──────────────────────────────────────────
    println("\n[1/8] Generating random sequences...")
    t_gen = elapsed() do
        return global BATCH = make_random_sequences(N_SEQUENCES, SEQ_LENGTH; seed=SEED)
    end
    println("      Done in $(fmt_ms(Float64(t_gen))) ms")
    println("      Batch: $(nsequences(BATCH)) sequences, length=$(seqlength(BATCH, 1)) bp")
    total_bp = N_SEQUENCES * SEQ_LENGTH
    println("      Total bases: $(total_bp)")

    # ── 2. Create models ──────────────────────────────────────────────────────
    println("\n[2/8] Creating models...")
    query_pwm = make_pwm(PWM_WIDTH, SEED)
    target_pwms = [make_pwm(PWM_WIDTH, SEED + 100 + i) for i in 1:N_TARGETS]
    println("      Query : $(query_pwm.name) (w=$PWM_WIDTH)")
    println("      Targets: $(N_TARGETS) distinct PWMs (w=$PWM_WIDTH)")

    # ── 3. Scan query model against sequences ──────────────────────────────────
    println("\n[3/8] Scanning query model (BestStrand, 10000 x 100)...")
    scan_times = bench_repeat(N_REPS) do
        return scan(query_pwm, BATCH; strands=BestStrand(), execution=execution)
    end
    s_scan = stats(scan_times)
    global QUERY_SCAN = scan(query_pwm, BATCH; strands=BestStrand(), execution=execution)
    println(
        "      min=$(fmt_ms(Float64(s_scan.min_ns))) ms  median=$(fmt_ms(Float64(s_scan.median_ns))) ms  max=$(fmt_ms(Float64(s_scan.max_ns))) ms",
    )

    # ── 4. Prepare query profile ───────────────────────────────────────────────
    println("\n[4/8] Preparing query profile (normalize + anchors)...")
    query_profile = ScoreProfile(query_pwm.name, QUERY_SCAN)
    prep_times = bench_repeat(N_REPS) do
        return prepare_profile(query_profile)
    end
    s_prep = stats(prep_times)
    global PREPARED = prepare_profile(query_profile)
    println(
        "      min=$(fmt_ms(Float64(s_prep.min_ns))) ms  median=$(fmt_ms(Float64(s_prep.median_ns))) ms  max=$(fmt_ms(Float64(s_prep.max_ns))) ms",
    )

    # ── 5. Scan one target model (timing) ──────────────────────────────────────
    println("\n[5/8] Scanning a single target model (BestStrand)...")
    single_scan_times = bench_repeat(N_REPS) do
        return scan(target_pwms[1], BATCH; strands=BestStrand(), execution=execution)
    end
    s_single = stats(single_scan_times)
    println(
        "      min=$(fmt_ms(Float64(s_single.min_ns))) ms  median=$(fmt_ms(Float64(s_single.median_ns))) ms  max=$(fmt_ms(Float64(s_single.max_ns))) ms",
    )

    # ── 6. Precompute all target profiles ──────────────────────────────────────
    println("\n[6/8] Scanning targets and constructing raw ScoreProfiles...")
    t_precompute = elapsed() do
        return global TARGET_PROFILES = [
            ScoreProfile(
                target_pwms[i].name,
                scan(target_pwms[i], BATCH; strands=BestStrand(), execution=execution),
            ) for i in 1:N_TARGETS
        ]
    end
    println(
        "      Done in $(fmt_ms(Float64(t_precompute))) ms ($(fmt_ms(Float64(t_precompute)/N_TARGETS)) ms/target)",
    )

    # ── 7. 1-vs-1 raw-profile comparison ──────────────────────────────────────
    println("\n[7/8] 1-vs-1 raw ScoreProfile comparison (normalize + align)...")
    one_times = bench_repeat(N_REPS) do
        return compare(
            PREPARED, TARGET_PROFILES[1]; metric=:co, search_range=10, window_radius=5
        )
    end
    s_one = stats(one_times)
    result1 = compare(
        PREPARED, TARGET_PROFILES[1]; metric=:co, search_range=10, window_radius=5
    )
    println(
        "      min=$(fmt_ms(Float64(s_one.min_ns))) ms  median=$(fmt_ms(Float64(s_one.median_ns))) ms  max=$(fmt_ms(Float64(s_one.max_ns))) ms",
    )
    println(
        "      Sample result: score=$(round(result1.score; digits=4)) offset=$(result1.offset) orient=$(result1.orientation) n_sites=$(result1.n_sites)",
    )

    # ── 8. 1-vs-50 raw-profile comparison ─────────────────────────────────────
    println("\n[8/8] 1-vs-50 raw ScoreProfile comparison (normalize + align)...")
    fifty_times = bench_repeat(N_REPS) do
        return compare(
            PREPARED,
            TARGET_PROFILES;
            execution=execution,
            metric=:co,
            search_range=10,
            window_radius=5,
        )
    end
    s_fifty = stats(fifty_times)
    results50 = compare(
        PREPARED,
        TARGET_PROFILES;
        execution=execution,
        metric=:co,
        search_range=10,
        window_radius=5,
    )
    println(
        "      min=$(fmt_ms(Float64(s_fifty.min_ns))) ms  median=$(fmt_ms(Float64(s_fifty.median_ns))) ms  max=$(fmt_ms(Float64(s_fifty.max_ns))) ms",
    )

    # ── End-to-end: scan + compare (query pre-prepared, targets scanned each time) ──
    println(
        "\n[E2E] Full end-to-end 1-vs-50 (query prepared, targets scanned + compared)..."
    )
    e2e_times = bench_repeat(N_REPS) do
        for t in target_pwms
            compare(
                PREPARED,
                t,
                BATCH;
                execution=execution,
                metric=:co,
                search_range=10,
                window_radius=5,
            )
        end
    end
    s_e2e = stats(e2e_times)
    println(
        "      min=$(fmt_ms(Float64(s_e2e.min_ns))) ms  median=$(fmt_ms(Float64(s_e2e.median_ns))) ms  max=$(fmt_ms(Float64(s_e2e.max_ns))) ms",
    )

    # ── Summary ────────────────────────────────────────────────────────────────
    println("\n" * "=" ^ 72)
    println("  SUMMARY")
    println("=" ^ 72)
    println(
        "  Sequence generation ($(N_SEQUENCES) x $(SEQ_LENGTH) bp):  $(fmt_ms(Float64(t_gen))) ms",
    )
    println(
        "  Query scan (10000 x 100, w=$PWM_WIDTH):        $(fmt_ms(Float64(s_scan.median_ns))) ms (median)",
    )
    println(
        "  Query profile preparation:              $(fmt_ms(Float64(s_prep.median_ns))) ms (median)",
    )
    println(
        "  Single target scan:                     $(fmt_ms(Float64(s_single.median_ns))) ms (median)",
    )
    println(
        "  Target scan + ScoreProfile construction: $(fmt_ms(Float64(t_precompute))) ms total",
    )
    println("  ---")
    println(
        "  1-vs-1  raw-profile comparison:          $(fmt_ms(Float64(s_one.median_ns))) ms (median)",
    )
    println(
        "  1-vs-50 raw-profile comparison:          $(fmt_ms(Float64(s_fifty.median_ns))) ms (median)",
    )
    println(
        "  End-to-end 1-vs-50 (scan+compare):       $(fmt_ms(Float64(s_e2e.median_ns))) ms (median)",
    )
    println("  ---")

    per_target_pure = Float64(s_fifty.median_ns) / N_TARGETS
    per_target_e2e = Float64(s_e2e.median_ns) / N_TARGETS
    println(
        "  Per-target (normalize + align):           $(@sprintf("%.3f", per_target_pure / 1e6)) ms/target",
    )
    println(
        "  Per-target (end-to-end):                 $(@sprintf("%.3f", per_target_e2e / 1e6)) ms/target",
    )
    println(
        "  Throughput (normalize + align):           $(@sprintf("%.0f", N_TARGETS / (Float64(s_fifty.median_ns) / 1e9))) comparisons/sec",
    )
    println(
        "  Throughput (end-to-end):                 $(@sprintf("%.0f", N_TARGETS / (Float64(s_e2e.median_ns) / 1e9))) comparisons/sec",
    )
    println(
        "  Batch speedup vs 50 scalar calls:         $(@sprintf("%.1f", (Float64(s_one.median_ns) * 50) / Float64(s_fifty.median_ns)))x",
    )

    # Show sample results
    println("\n  Sample 1-vs-50 results (first 5):")
    for i in 1:min(5, length(results50))
        r = results50[i]
        println(
            "    $(r.target): score=$(@sprintf("%.4f", r.score)) offset=$(r.offset) orient=$(r.orientation) n_sites=$(r.n_sites)",
        )
    end
    println("=" ^ 72)

    return nothing
end

main()
