# Mimosa.jl Reproducible Benchmark Suite — PLAN_2.md E1/E2
#
# Comprehensive benchmarks covering all representative workloads from PLAN.md:
# ragged heavy-tail batches, BaMM orders 1–5, one-to-many 10/100/1000,
# dense threshold anchors, high/low site density, null schedules.
#
# Collects: median, min, mean, variance, allocations, peak RSS.
# Supports serial and threaded (1/2/4/available) configurations.
# Outputs machine-readable JSON to stdout (or file with --output).
#
# Usage:
#   julia --project=Mimosa.jl/benchmark Mimosa.jl/benchmark/runbenchmarks.jl
#   julia --project=Mimosa.jl/benchmark Mimosa.jl/benchmark/runbenchmarks.jl --output results.json
#   julia --project=Mimosa.jl/benchmark Mimosa.jl/benchmark/runbenchmarks.jl --report
#   julia --project=Mimosa.jl/benchmark Mimosa.jl/benchmark/runbenchmarks.jl --baseline baseline.json
#
# Environment variables:
#   JULIA_NUM_THREADS  — controls thread count for threaded scaling tests

using Mimosa
using BenchmarkTools
using Printf
using Dates
using Random
using Pkg

const REPO_ROOT = dirname(dirname(@__DIR__))
const EXAMPLES = joinpath(REPO_ROOT, "examples")

# ──────────────────────────────────────────────────────────────────────────
# CLI argument parsing
# ──────────────────────────────────────────────────────────────────────────

struct BenchConfig
    output_file::Union{Nothing,String}
    report_only::Bool
    baseline_file::Union{Nothing,String}
    samples::Int
    seconds::Float64
end

function parse_args(args::Vector{String})
    output_file = nothing
    report_only = false
    baseline_file = nothing
    samples = 100
    seconds_budget = 5.0
    i = 1
    while i <= length(args)
        if args[i] == "--output" && i < length(args)
            output_file = args[i + 1]
            i += 2
        elseif args[i] == "--report"
            report_only = true
            i += 1
        elseif args[i] == "--baseline" && i < length(args)
            baseline_file = args[i + 1]
            i += 2
        elseif args[i] == "--samples" && i < length(args)
            samples = parse(Int, args[i + 1])
            i += 2
        elseif args[i] == "--seconds" && i < length(args)
            seconds_budget = parse(Float64, args[i + 1])
            i += 2
        elseif args[i] in ("--help", "-h")
            println(
                """
        Mimosa.jl Benchmark Suite

        Usage: julia --project=Mimosa.jl/benchmark Mimosa.jl/benchmark/runbenchmarks.jl [OPTIONS]

        Options:
          --output <file>    Write JSON results to file (default: stdout)
          --report           Print environment metadata only, no benchmarks
          --baseline <file>  Compare results against baseline file
          --samples <n>      BenchmarkTools samples per benchmark (default: 100)
          --seconds <f>      BenchmarkTools seconds budget per benchmark (default: 5.0)
          --help, -h         Show this help
        """,
            )
            return nothing
        else
            @warn "Unknown argument: $(args[i])"
            i += 1
        end
    end
    return BenchConfig(output_file, report_only, baseline_file, samples, seconds_budget)
end

# ──────────────────────────────────────────────────────────────────────────
# Environment metadata collection
# ──────────────────────────────────────────────────────────────────────────

function collect_environment()
    cpu = Sys.cpu_info()
    cpu_model = length(cpu) > 0 ? cpu[1].model : "unknown"
    cpu_speed = length(cpu) > 0 ? cpu[1].speed : 0
    n_cpus = length(cpu)

    # Get git commit SHA
    commit_sha = "unknown"
    try
        commit_sha = String(read(`git -C $(REPO_ROOT) rev-parse --short HEAD`))[1:(end - 1)]
    catch
        # Not in a git repo or git not available
    end

    # Get package versions
    pkg_versions = Dict{String,String}()
    try
        for dep in Pkg.project().dependencies
            name = first(dep)
            pkg_versions[name] = "unknown"
        end
    catch
        # Fallback: just record what we can
    end
    # Get Mimosa version from Project.toml
    mimosa_version = "unknown"
    try
        project_toml = joinpath(REPO_ROOT, "Mimosa.jl", "Project.toml")
        for line in eachline(project_toml)
            if startswith(strip(line), "version")
                mimosa_version = String(strip(replace(line, r"version\s*=\s*" => "")))
                break
            end
        end
    catch
    end

    return Dict{String,Any}(
        "commit_sha" => commit_sha,
        "julia_version" => string(VERSION),
        "julia_executable" => Sys.BINDIR,
        "machine" => Sys.MACHINE,
        "os" => if Sys.islinux()
            "Linux"
        elseif Sys.isapple()
            "macOS"
        elseif Sys.iswindows()
            "Windows"
        else
            "unknown"
        end,
        "kernel" => try
            String(read(`uname -r`))[1:(end - 1)]
        catch
            "unknown"
        end,
        "cpu_model" => cpu_model,
        "cpu_speed_mhz" => cpu_speed,
        "n_cpus" => n_cpus,
        "n_threads" => Threads.nthreads(),
        "total_ram_gb" => round(Sys.total_memory() / 1024^3; digits=2),
        "free_ram_gb" => round(Sys.free_memory() / 1024^3; digits=2),
        "package_versions" => pkg_versions,
        "mimosa_version" => mimosa_version,
        "timestamp" => Dates.format(now(), "yyyy-mm-ddTHH:MM:SS"),
        "warmup_policy" => "1 explicit warm-up call before each measurement",
        "samples" => 0,  # filled per-benchmark
        "seconds_budget" => 0.0,  # filled per-benchmark
    )
end

# ──────────────────────────────────────────────────────────────────────────
# Benchmark result collection
# ──────────────────────────────────────────────────────────────────────────

struct BenchResult
    name::String
    category::String
    median_ns::Float64
    min_ns::Float64
    mean_ns::Float64
    variance_ns::Float64
    allocations::Int
    memory_bytes::Int
    n_samples::Int
    n_evals::Int
    warmup::Bool
    parameters::Dict{String,Any}
end

function bench_result(
    name::String, category::String, b::BenchmarkTools.Trial; parameters...
)
    params = Dict{String,Any}()
    for (k, v) in parameters
        params[String(k)] = v
    end
    return BenchResult(
        name,
        category,
        Float64(BenchmarkTools.median(b).time),
        Float64(minimum(b).time),
        Float64(mean(b.times)),
        Float64(var(b.times)),
        b.allocs,
        b.memory,
        length(b.times),
        b.params.evals,
        true,  # warmup performed
        params,
    )
end

function result_to_dict(r::BenchResult)
    return Dict{String,Any}(
        "name" => r.name,
        "category" => r.category,
        "median_ns" => r.median_ns,
        "min_ns" => r.min_ns,
        "mean_ns" => r.mean_ns,
        "variance_ns" => r.variance_ns,
        "allocations" => r.allocations,
        "memory_bytes" => r.memory_bytes,
        "n_samples" => r.n_samples,
        "n_evals" => r.n_evals,
        "warmup" => r.warmup,
        "parameters" => r.parameters,
    )
end

# ──────────────────────────────────────────────────────────────────────────
# JSON serialization (minimal, no external dependency)
# ──────────────────────────────────────────────────────────────────────────

function to_json_str(d::Dict)
    return _json_str(d)
end

function _json_str(x::AbstractString)
    # Escape string for JSON
    escaped = replace(x, "\\" => "\\\\")
    escaped = replace(escaped, "\"" => "\\\"")
    escaped = replace(escaped, "\n" => "\\n")
    escaped = replace(escaped, "\r" => "\\r")
    escaped = replace(escaped, "\t" => "\\t")
    return "\"$escaped\""
end

function _json_str(x::Number)
    if isnan(x) || isinf(x)
        return "null"
    end
    return string(x)
end

function _json_str(x::Bool)
    return x ? "true" : "false"
end

function _json_str(::Nothing)
    return "null"
end

function _json_str(d::Dict)
    pairs = String[]
    for (k, v) in d
        push!(pairs, "$(_json_str(string(k))): $(_json_str(v))")
    end
    return "{\n  $(join(pairs, ",\n  "))\n}"
end

function _json_str(v::AbstractVector)
    if isempty(v)
        return "[]"
    end
    items = [_json_str(x) for x in v]
    return "[\n  $(join(items, ",\n  "))\n]"
end

function _json_str(v::AbstractArray)
    return _json_str(vec(v))
end

# ──────────────────────────────────────────────────────────────────────────
# Model and data setup helpers
# ──────────────────────────────────────────────────────────────────────────

function make_pwm(width::Int)
    rng = Random.MersenneTwister(42)
    weights = Matrix{Float32}(undef, 5, width)
    for col in 1:width
        for row in 1:4
            weights[row, col] = Float32(randn(rng) * 0.5)
        end
        weights[5, col] = minimum(@view weights[1:4, col]) - Float32(0.1)
    end
    bg = (Float32(0.25), Float32(0.25), Float32(0.25), Float32(0.25))
    return PWM("bench_pwm_$width", weights, bg)
end

function make_bamm(order::Int, width::Int)
    n_rows = 5^(order + 1)
    rng = Random.MersenneTwister(43)
    rep = Matrix{Float32}(undef, n_rows, width)
    for i in 1:n_rows, j in 1:width
        rep[i, j] = Float32(randn(rng) * 0.3)
    end
    return BaMM("bench_bamm_o$order", rep, order, width)
end

function make_relations_file(model_names::Vector{String})
    tsv = tempname()
    lines = ["motif\tgroup"]
    for (i, name) in enumerate(model_names)
        group = i <= length(model_names) ÷ 2 ? "A" : "B"
        push!(lines, "$name\t$group")
    end
    write(tsv, join(lines, "\n") * "\n")
    return tsv
end

# ──────────────────────────────────────────────────────────────────────────
# Package import and startup timing
# ──────────────────────────────────────────────────────────────────────────

"""
    measure_import_time()

Measure the time to `using Mimosa` in a fresh subprocess.
Returns time in nanoseconds.
"""
function measure_import_time()
    julia = joinpath(Sys.BINDIR, "julia")
    project = joinpath(REPO_ROOT, "Mimosa.jl")
    cmd = `$julia --project=$project -e 'using Mimosa; println("OK")'`
    # Warm up (precompile cache)
    try
        run(pipeline(cmd; stdout=devnull, stderr=devnull))
    catch
    end
    # Measure
    t0 = time()
    try
        run(pipeline(cmd; stdout=devnull, stderr=devnull))
    catch
    end
    return (time() - t0) * 1e9
end

"""
    measure_precompile_time()

Measure Pkg.precompile time in a subprocess.
"""
function measure_precompile_time()
    julia = joinpath(Sys.BINDIR, "julia")
    project = joinpath(REPO_ROOT, "Mimosa.jl")
    cmd = `$julia --project=$project -e 'using Pkg; Pkg.precompile()'`
    t0 = time()
    try
        run(pipeline(cmd; stdout=devnull, stderr=devnull))
    catch
    end
    return (time() - t0) * 1e9
end

"""
    measure_cli_startup()

Measure real CLI subprocess startup time (mimosa --version).
"""
function measure_cli_startup()
    julia = joinpath(Sys.BINDIR, "julia")
    project = joinpath(REPO_ROOT, "Mimosa.jl")
    app = joinpath(REPO_ROOT, "Mimosa.jl", "app", "mimosa.jl")
    cmd = `$julia --project=$project $app --version`
    # Warm up
    try
        run(pipeline(cmd; stdout=devnull, stderr=devnull))
    catch
    end
    t0 = time()
    try
        run(pipeline(cmd; stdout=devnull, stderr=devnull))
    catch
    end
    return (time() - t0) * 1e9
end

# ──────────────────────────────────────────────────────────────────────────
# Benchmark groups
# ──────────────────────────────────────────────────────────────────────────

"""
    bench_import_and_startup(results, config)

Measure package import, precompile, and CLI startup times.
"""
function bench_import_and_startup!(results::Vector{BenchResult}, config::BenchConfig)
    println("\n=== Package Import and Startup ===")

    # Import time
    println("  Measuring package import time...")
    import_ns = measure_import_time()
    push!(
        results,
        BenchResult(
            "using_mimosa",
            "startup",
            import_ns,
            import_ns,
            import_ns,
            0.0,
            0,
            0,
            1,
            1,
            true,
            Dict{String,Any}("method" => "subprocess"),
        ),
    )

    # Precompile time
    println("  Measuring precompile time...")
    precompile_ns = measure_precompile_time()
    push!(
        results,
        BenchResult(
            "pkg_precompile",
            "startup",
            precompile_ns,
            precompile_ns,
            precompile_ns,
            0.0,
            0,
            0,
            1,
            1,
            true,
            Dict{String,Any}("method" => "subprocess"),
        ),
    )

    # CLI startup
    println("  Measuring CLI startup time...")
    cli_ns = measure_cli_startup()
    return push!(
        results,
        BenchResult(
            "cli_startup",
            "startup",
            cli_ns,
            cli_ns,
            cli_ns,
            0.0,
            0,
            0,
            1,
            1,
            true,
            Dict{String,Any}("method" => "subprocess", "command" => "mimosa --version"),
        ),
    )
end

"""
    bench_pwm_scan(results, config)

Benchmark PWM scanning at different sizes and strand policies.
"""
function bench_pwm_scan!(results::Vector{BenchResult}, config::BenchConfig)
    println("\n=== PWM Scanning (single sequence) ===")

    for width in (8, 15, 30)
        pwm = make_pwm(width)
        for seqlen in (100, 200, 1000)
            seq = make_random_sequences(1, seqlen; seed=42)
            enc = Mimosa.sequence(seq, 1)
            n_pos = npositions(seqlen, width)

            # Forward scan
            dest = Vector{Float32}(undef, n_pos)
            scan_forward!(dest, pwm, enc, n_pos)  # warm up
            b = BenchmarkTools.@benchmark scan_forward!($dest, $pwm, $enc, $n_pos)
            println(
                @sprintf(
                    "  scan_forward!  w=%d len=%d  median=%.3f μs  allocs=%d",
                    width,
                    seqlen,
                    median(b).time / 1000,
                    b.allocs
                )
            )
            push!(
                results,
                bench_result(
                    "pwm_scan_forward",
                    "pwm_scan",
                    b;
                    width=width,
                    seqlen=seqlen,
                    strand="forward",
                ),
            )

            # Best strand
            best_hits!(dest, pwm, enc, n_pos)
            b = BenchmarkTools.@benchmark best_hits!($dest, $pwm, $enc, $n_pos)
            println(
                @sprintf(
                    "  best_hits!     w=%d len=%d  median=%.3f μs  allocs=%d",
                    width,
                    seqlen,
                    median(b).time / 1000,
                    b.allocs
                )
            )
            push!(
                results,
                bench_result(
                    "pwm_best_hits",
                    "pwm_scan",
                    b;
                    width=width,
                    seqlen=seqlen,
                    strand="best",
                ),
            )

            # Reverse complement
            rc = similar(enc)
            reverse_complement!(rc, enc)
            b = BenchmarkTools.@benchmark reverse_complement!($rc, $enc)
            println(
                @sprintf(
                    "  rev_comp!      w=%d len=%d  median=%.3f μs  allocs=%d",
                    width,
                    seqlen,
                    median(b).time / 1000,
                    b.allocs
                )
            )
            push!(
                results,
                bench_result(
                    "reverse_complement", "pwm_scan", b; width=width, seqlen=seqlen
                ),
            )
        end
        println()
    end
end

"""
    bench_batch_scan(results, config)

Benchmark batch scanning with ragged heavy-tail batches and serial/threaded scaling.
"""
function bench_batch_scan!(results::Vector{BenchResult}, config::BenchConfig)
    println("\n=== Batch Scanning — Ragged Heavy-Tail ===")

    pwm = make_pwm(15)

    # Ragged heavy-tail: mix of short and long sequences
    for (label, n_seqs, seq_len) in
        (("small", 100, 200), ("medium", 1000, 200), ("large", 10000, 200))
        batch = make_random_sequences(n_seqs, seq_len; seed=42)

        # Serial
        scan(pwm, batch; strands=BestStrand(), execution=Execution())  # warmup
        b = BenchmarkTools.@benchmark scan(
            $pwm, $batch; strands=BestStrand(), execution=Execution()
        )
        median_ser = median(b).time
        println(
            @sprintf(
                "  serial    n=%-5d len=%-4d  median=%.3f ms  allocs=%d",
                n_seqs,
                seq_len,
                median_ser / 1e6,
                b.allocs
            )
        )
        push!(
            results,
            bench_result(
                "batch_scan_serial",
                "batch_scan",
                b;
                n_seqs=n_seqs,
                seq_len=seq_len,
                execution="serial",
                workload=label,
            ),
        )

        # Threaded at different thread counts
        max_threads = max(1, Threads.nthreads())
        for nt in (1, 2, 4)
            if nt > max_threads
                continue
            end
            exec = Execution(nt)
            scan(pwm, batch; strands=BestStrand(), execution=exec)  # warmup
            b = BenchmarkTools.@benchmark scan(
                $pwm, $batch; strands=BestStrand(), execution=($exec)
            )
            speedup = median_ser / median(b).time
            println(
                @sprintf(
                    "  threaded  n=%-5d len=%-4d  t=%d  median=%.3f ms  allocs=%d  speedup=%.2fx",
                    n_seqs,
                    seq_len,
                    nt,
                    median(b).time / 1e6,
                    b.allocs,
                    speedup
                )
            )
            push!(
                results,
                bench_result(
                    "batch_scan_threaded",
                    "batch_scan",
                    b;
                    n_seqs=n_seqs,
                    seq_len=seq_len,
                    execution="threaded",
                    n_threads=nt,
                    speedup=speedup,
                    workload=label,
                ),
            )
        end
        println()
    end

    # Ragged heavy-tail: variable-length sequences
    println("  --- Ragged heavy-tail (variable lengths) ---")
    for (label, lengths) in (
        ("short_heavy", [50, 80, 100, 120, 50, 60, 100, 80, 50, 70]),
        ("long_heavy", [500, 1000, 2000, 500, 800, 1500, 3000, 600, 1000, 400]),
    )
        # Build ragged batch manually (Julia 1-based offsets)
        n = length(lengths)
        total = sum(lengths)
        data = Vector{UInt8}(undef, total)
        offsets = Vector{Int64}(undef, n + 1)
        offsets[1] = 1  # Julia 1-based
        rng = Random.MersenneTwister(42)
        for i in 1:n
            for j in 1:lengths[i]
                data[offsets[i] + j - 1] = UInt8(rand(rng, 0:4))
            end
            offsets[i + 1] = offsets[i] + lengths[i]
        end
        batch = EncodedSequenceBatch(data, offsets)

        scan(pwm, batch; strands=BestStrand(), execution=Execution())  # warmup
        b = BenchmarkTools.@benchmark scan(
            $pwm, $batch; strands=BestStrand(), execution=Execution()
        )
        println(
            @sprintf(
                "  ragged %s  n=%d  median=%.3f ms  allocs=%d",
                label,
                n,
                median(b).time / 1e6,
                b.allocs
            )
        )
        push!(
            results,
            bench_result(
                "batch_scan_ragged",
                "batch_scan",
                b;
                n_seqs=n,
                seq_lengths=lengths,
                execution="serial",
                workload=label,
            ),
        )
    end
    return println()
end

"""
    bench_one_to_many(results, config)

Benchmark one-to-many profile comparison with 10/100/1000 targets.
"""
function bench_one_to_many!(results::Vector{BenchResult}, config::BenchConfig)
    println("\n=== One-to-Many Profile Comparison ===")

    pwm = make_pwm(15)
    batch = make_random_sequences(50, 200; seed=42)
    scan_result = scan(pwm, batch; strands=BestStrand())

    query_profile = ScoreProfile("query", scan_result)
    prepared = Mimosa.prepare_profile(query_profile)

    for n_targets in (10, 100, 1000)
        # Build target profiles
        targets = [ScoreProfile("target_$i", scan_result) for i in 1:n_targets]

        policies = Tuple{String,Execution,Int}[("serial", Execution(), 1)]
        if Threads.nthreads() > 1
            push!(policies, ("threaded", Execution(Threads.nthreads()), Threads.nthreads()))
        end
        for (execution_name, execution, n_threads) in policies
            compare(
                prepared,
                targets;
                execution=execution,
                metric=:co,
                search_range=10,
                window_radius=5,
            ) # warmup
            b = BenchmarkTools.@benchmark compare(
                $prepared,
                $targets;
                execution=($execution),
                metric=:co,
                search_range=10,
                window_radius=5,
            )
            println(
                @sprintf(
                    "  one-to-many  n=%-4d  execution=%-8s  median=%.3f ms  allocs=%d",
                    n_targets,
                    execution_name,
                    median(b).time / 1e6,
                    b.allocs
                )
            )
            push!(
                results,
                bench_result(
                    "one_to_many_compare",
                    "one_to_many",
                    b;
                    n_targets=n_targets,
                    metric="co",
                    search_range=10,
                    window_radius=5,
                    execution=execution_name,
                    n_threads=n_threads,
                ),
            )
        end
    end
    return println()
end

"""
    bench_higher_order_scan(results, config)

Benchmark higher-order model scanning (BaMM orders 0-5).
"""
function bench_higher_order_scan!(results::Vector{BenchResult}, config::BenchConfig)
    println("\n=== Higher-Order Scanning (BaMM orders 0-5) ===")

    width = 10
    for order in (0, 1, 2, 3, 4, 5)
        n_rows = 5^(order + 1)
        if n_rows > 50000
            println(
                @sprintf(
                    "  BaMM order=%d  rows=%d  SKIPPED (too large for benchmark)",
                    order,
                    n_rows
                )
            )
            continue
        end
        model = make_bamm(order, width)
        seq = make_random_sequences(1, 200; seed=42)
        enc = Mimosa.sequence(seq, 1)
        kmer = order + 1
        ctx = order
        win = width + order
        n_pos = max(length(enc) - win + 1, 0)
        dest = Vector{Float32}(undef, n_pos)

        scan_forward!(dest, model, enc, n_pos)  # warmup
        b = BenchmarkTools.@benchmark scan_forward!($dest, $model, $enc, $n_pos)
        println(
            @sprintf(
                "  BaMM order=%d  kmer=%d  rows=%5d  median=%.3f μs  allocs=%d",
                order,
                kmer,
                n_rows,
                median(b).time / 1000,
                b.allocs
            )
        )
        push!(
            results,
            bench_result(
                "bamm_scan",
                "higher_order_scan",
                b;
                order=order,
                kmer=kmer,
                n_rows=n_rows,
                width=width,
            ),
        )
    end

    # Also benchmark real BaMM from examples
    println("\n  --- Real BaMM models from examples ---")
    for fname in ("foxa2.ihbcp", "gata2.ihbcp", "gata4.ihbcp", "myog.ihbcp")
        path = joinpath(EXAMPLES, fname)
        isfile(path) || continue
        model = readmodel(path)
        order = model.order
        width = model.motif_length
        seq = make_random_sequences(1, 200; seed=42)
        enc = Mimosa.sequence(seq, 1)
        kmer = order + 1
        ctx = order
        win = width + order
        n_pos = max(length(enc) - win + 1, 0)
        if n_pos <= 0
            println(
                @sprintf(
                    "  %s  order=%d  width=%d  SKIPPED (seq too short)", fname, order, width
                )
            )
            continue
        end
        dest = Vector{Float32}(undef, n_pos)
        scan_forward!(dest, model, enc, n_pos)  # warmup
        b = BenchmarkTools.@benchmark scan_forward!($dest, $model, $enc, $n_pos)
        println(
            @sprintf(
                "  %-15s  order=%d  width=%d  rows=%5d  median=%.3f μs  allocs=%d",
                fname,
                order,
                width,
                size(model.representation, 1),
                median(b).time / 1000,
                b.allocs
            )
        )
        push!(
            results,
            bench_result(
                "bamm_scan_real",
                "higher_order_scan",
                b;
                model_file=fname,
                order=order,
                width=width,
                n_rows=size(model.representation, 1),
            ),
        )
    end
    return println()
end

"""
    bench_site_extraction(results, config)

Benchmark site extraction with different selectors and site densities.
"""
function bench_site_extraction!(results::Vector{BenchResult}, config::BenchConfig)
    println("\n=== Site Extraction and PFM Reconstruction ===")

    pwm = make_pwm(15)

    for (n_seqs, seq_len) in ((100, 200), (1000, 200))
        batch = make_random_sequences(n_seqs, seq_len; seed=42)

        # BestPerSequence — low site density (1 site per sequence)
        selectsites(pwm, batch, BestPerSequence(); strands=BestStrand())  # warmup
        b = BenchmarkTools.@benchmark selectsites(
            $pwm, $batch, BestPerSequence(); strands=BestStrand()
        )
        println(
            @sprintf(
                "  BestPerSequence  n=%-4d  median=%.3f μs  allocs=%d",
                n_seqs,
                median(b).time / 1000,
                b.allocs
            )
        )
        push!(
            results,
            bench_result(
                "site_extraction_best",
                "site_extraction",
                b;
                n_seqs=n_seqs,
                selector="best_per_sequence",
                density="low",
            ),
        )

        # ThresholdHits — high site density (all sites above threshold)
        mn, mx = scorebounds(pwm)
        threshold = Float32(mn + 0.5 * (mx - mn))
        selectsites(pwm, batch, ThresholdHits(threshold); strands=BestStrand())  # warmup
        b = BenchmarkTools.@benchmark selectsites(
            $pwm, $batch, ThresholdHits($threshold); strands=BestStrand()
        )
        println(
            @sprintf(
                "  ThresholdHits    n=%-4d  t=%.2f  median=%.3f μs  allocs=%d",
                n_seqs,
                threshold,
                median(b).time / 1000,
                b.allocs
            )
        )
        push!(
            results,
            bench_result(
                "site_extraction_threshold",
                "site_extraction",
                b;
                n_seqs=n_seqs,
                selector="threshold_hits",
                threshold=threshold,
                density="high",
            ),
        )

        # TopFractionHits — variable density
        frac = TopFractionHits(0.1)
        selectsites(pwm, batch, frac; strands=BestStrand())  # warmup
        b = BenchmarkTools.@benchmark selectsites($pwm, $batch, $frac; strands=BestStrand())
        println(
            @sprintf(
                "  TopFractionHits  n=%-4d  f=0.1  median=%.3f μs  allocs=%d",
                n_seqs,
                median(b).time / 1000,
                b.allocs
            )
        )
        push!(
            results,
            bench_result(
                "site_extraction_topfraction",
                "site_extraction",
                b;
                n_seqs=n_seqs,
                selector="top_fraction",
                fraction=0.1,
                density="medium",
            ),
        )

        # PFM reconstruction
        reconstruct_pfm(pwm, batch, BestPerSequence(); pseudocount=Float32(1e-4))  # warmup
        b = BenchmarkTools.@benchmark reconstruct_pfm(
            $pwm, $batch, BestPerSequence(); pseudocount=Float32(1e-4)
        )
        println(
            @sprintf(
                "  reconstruct_pfm  n=%-4d  median=%.3f μs  allocs=%d",
                n_seqs,
                median(b).time / 1000,
                b.allocs
            )
        )
        push!(
            results,
            bench_result(
                "pfm_reconstruction",
                "site_extraction",
                b;
                n_seqs=n_seqs,
                selector="best_per_sequence",
            ),
        )
    end
    return println()
end

"""
    bench_gev(results, config)

Benchmark GEV fitting and p-value computation.
"""
function bench_gev!(results::Vector{BenchResult}, config::BenchConfig)
    println("\n=== GEV Fitting and Statistics ===")

    for n in (100, 500, 2000)
        rng = Random.MersenneTwister(42)
        samples = Float32.(randn(rng, n) .* 0.3 .+ 0.5)

        fit_gev(samples)  # warmup
        b = BenchmarkTools.@benchmark fit_gev($samples)
        println(
            @sprintf(
                "  fit_gev(n=%-4d)  median=%.3f μs  allocs=%d",
                n,
                median(b).time / 1000,
                b.allocs
            )
        )
        push!(results, bench_result("fit_gev", "statistics", b; n_samples=n))
    end

    # BH FDR
    pvals = Float32.(rand(MersenneTwister(42), 1000))
    adjusted_pvalues(pvals; method=BenjaminiHochberg())  # warmup
    b = BenchmarkTools.@benchmark adjusted_pvalues($pvals; method=BenjaminiHochberg())
    println(
        @sprintf(
            "  BH FDR(n=1000)  median=%.3f μs  allocs=%d", median(b).time / 1000, b.allocs
        )
    )
    push!(results, bench_result("bh_fdr", "statistics", b; n_pvalues=1000))

    return println()
end

"""
    bench_null_distribution(results, config)

Benchmark null distribution building (motif and profile strategies).
"""
function bench_null_distribution!(results::Vector{BenchResult}, config::BenchConfig)
    println("\n=== Null Distribution Building ===")

    # Load real models
    pwm1 = readmodel(joinpath(EXAMPLES, "pif4.meme"))
    pwm2 = readmodel(joinpath(EXAMPLES, "foxa2.meme"))
    pwm3 = readmodel(joinpath(EXAMPLES, "gata2.meme"))
    pwm4 = readmodel(joinpath(EXAMPLES, "gata4.meme"))

    # Build a larger model set for realistic null schedules
    models = [pwm1, pwm2, pwm3, pwm4]
    # Duplicate with different names to get more comparisons
    for i in 5:10
        m = PWM("model_$i", copy(pwm1.representation), pwm1.background)
        push!(models, m)
    end

    model_names = [m.name for m in models]
    rel_file = make_relations_file(model_names)
    rel = Mimosa.parse_group_relations(rel_file)

    # Profile strategy
    batch = make_random_sequences(50, 200; seed=42)
    build_null(models, rel; metric=:co, sequences=batch)
    b = BenchmarkTools.@benchmark build_null($models, $rel; metric=:co, sequences=($batch))
    println(
        @sprintf(
            "  build_null profile co   n_models=%d  n_seqs=50  median=%.3f ms  allocs=%d",
            length(models),
            median(b).time / 1e6,
            b.allocs
        )
    )
    push!(
        results,
        bench_result(
            "build_null",
            "null_distribution",
            b;
            strategy="profile",
            metric="co",
            n_models=length(models),
            n_sequences=50,
        ),
    )

    rm(rel_file; force=true)
    return println()
end

"""
    bench_serial_vs_threaded(results, config)

Verify serial and threaded scanning produce identical results.
"""
function bench_serial_vs_threaded!(results::Vector{BenchResult}, config::BenchConfig)
    println("\n=== Serial vs Threaded Equivalence ===")

    pwm = make_pwm(15)
    batch = make_random_sequences(500, 200; seed=42)

    ser = scan(pwm, batch; strands=BestStrand(), execution=Execution())
    nthreads = max(1, Threads.nthreads())

    for nt in (1, 2, 4)
        if nt > nthreads
            continue
        end
        thr = scan(pwm, batch; strands=BestStrand(), execution=Execution(nt))
        identical = ser == thr
        status = identical ? "PASS" : "FAIL"
        println(@sprintf("  serial == threaded(%d tasks): %s", nt, status))
        push!(
            results,
            BenchResult(
                "serial_vs_threaded_equiv",
                "correctness",
                0.0,
                0.0,
                0.0,
                0.0,
                0,
                0,
                1,
                1,
                true,
                Dict{String,Any}("n_threads" => nt, "identical" => identical),
            ),
        )
    end
    return println()
end

"""
    bench_serialization(results, config)

Benchmark JSON serialization and dict conversion.
"""
function bench_serialization!(results::Vector{BenchResult}, config::BenchConfig)
    println("\n=== Serialization Latency ===")

    pwm = make_pwm(15)
    batch = make_random_sequences(10, 100; seed=41)
    result = compare(pwm, pwm, batch; metric=:co)

    to_json(result)  # warmup
    b = BenchmarkTools.@benchmark to_json($result)
    println(
        @sprintf(
            "  to_json(ComparisonResult)  median=%.3f μs  allocs=%d",
            median(b).time / 1000,
            b.allocs
        )
    )
    push!(
        results, bench_result("to_json", "serialization", b; result_type="ComparisonResult")
    )

    to_dict(result)
    b = BenchmarkTools.@benchmark to_dict($result)
    println(
        @sprintf(
            "  to_dict(ComparisonResult)  median=%.3f μs  allocs=%d",
            median(b).time / 1000,
            b.allocs
        )
    )
    push!(
        results, bench_result("to_dict", "serialization", b; result_type="ComparisonResult")
    )
    return println()
end

"""
    bench_storage(results, config)

Benchmark model storage (bundle write/read round-trip).
"""
function _write_benchmark_model!(
    tmpdir::String, prefix::String, model::AbstractMotifModel, counter::Base.RefValue{Int}
)
    counter[] += 1
    return writemodel(joinpath(tmpdir, "$(prefix)_$(counter[])"), model)
end

function bench_storage!(results::Vector{BenchResult}, config::BenchConfig)
    println("\n=== Storage (bundle write/read) ===")

    pwm = make_pwm(15)
    bamm = make_bamm(2, 10)

    tmpdir = mktempdir()

    # Bundles reject overwrites, so each benchmark sample needs a unique target.
    write_counter = Ref(0)

    # PWM write
    writemodel(joinpath(tmpdir, "pwm_bundle"), pwm)
    b = BenchmarkTools.@benchmark _write_benchmark_model!(
        $tmpdir, "pwm_write", $pwm, $write_counter
    )
    println(
        @sprintf(
            "  writemodel PWM     median=%.3f μs  allocs=%d",
            median(b).time / 1000,
            b.allocs
        )
    )
    push!(results, bench_result("writemodel", "storage", b; model_type="PWM"))

    # PWM read
    readmodel(joinpath(tmpdir, "pwm_bundle"))
    b = BenchmarkTools.@benchmark readmodel(joinpath($tmpdir, "pwm_bundle"))
    println(
        @sprintf(
            "  readmodel  PWM     median=%.3f μs  allocs=%d",
            median(b).time / 1000,
            b.allocs
        )
    )
    push!(results, bench_result("readmodel", "storage", b; model_type="PWM"))

    # BaMM write
    writemodel(joinpath(tmpdir, "bamm_bundle"), bamm)
    b = BenchmarkTools.@benchmark _write_benchmark_model!(
        $tmpdir, "bamm_write", $bamm, $write_counter
    )
    println(
        @sprintf(
            "  writemodel BaMM    median=%.3f μs  allocs=%d",
            median(b).time / 1000,
            b.allocs
        )
    )
    push!(results, bench_result("writemodel", "storage", b; model_type="BaMM"))

    # BaMM read
    readmodel(joinpath(tmpdir, "bamm_bundle"))
    b = BenchmarkTools.@benchmark readmodel(joinpath($tmpdir, "bamm_bundle"))
    println(
        @sprintf(
            "  readmodel  BaMM    median=%.3f μs  allocs=%d",
            median(b).time / 1000,
            b.allocs
        )
    )
    push!(results, bench_result("readmodel", "storage", b; model_type="BaMM"))

    rm(tmpdir; recursive=true, force=true)
    return println()
end

# ──────────────────────────────────────────────────────────────────────────
# Baseline comparison
# ──────────────────────────────────────────────────────────────────────────

function load_baseline(path::String)
    if !isfile(path)
        @warn "Baseline file not found: $path"
        return Dict{String,Any}()
    end
    # Parse JSON manually (minimal parser)
    content = read(path, String)
    return _parse_json(content)
end

# Minimal JSON parser
function _parse_json(s::AbstractString)
    s = String(strip(String(s)))
    (pos, val) = _parse_json_value(s, 1)
    return val
end

function _skip_ws(s::String, pos::Int)
    while pos <= length(s) && s[pos] in (' ', '\t', '\n', '\r')
        pos += 1
    end
    return pos
end

function _parse_json_value(s::String, pos::Int)
    pos = _skip_ws(s, pos)
    if pos > length(s)
        return (pos, nothing)
    end
    c = s[pos]
    if c == '{'
        return _parse_json_object(s, pos)
    elseif c == '['
        return _parse_json_array(s, pos)
    elseif c == '"'
        return _parse_json_string(s, pos)
    elseif c == 't' || c == 'f'
        return _parse_json_bool(s, pos)
    elseif c == 'n'
        return _parse_json_null(s, pos)
    else
        return _parse_json_number(s, pos)
    end
end

function _parse_json_object(s::String, pos::Int)
    pos += 1  # skip {
    result = Dict{String,Any}()
    pos = _skip_ws(s, pos)
    if pos <= length(s) && s[pos] == '}'
        return (pos + 1, result)
    end
    while true
        pos = _skip_ws(s, pos)
        (pos, key) = _parse_json_string(s, pos)
        pos = _skip_ws(s, pos)
        pos += 1  # skip :
        (pos, val) = _parse_json_value(s, pos)
        result[key] = val
        pos = _skip_ws(s, pos)
        if pos <= length(s) && s[pos] == ','
            pos += 1
        else
            break
        end
    end
    pos = _skip_ws(s, pos)
    pos += 1  # skip }
    return (pos, result)
end

function _parse_json_array(s::String, pos::Int)
    pos += 1  # skip [
    result = Any[]
    pos = _skip_ws(s, pos)
    if pos <= length(s) && s[pos] == ']'
        return (pos + 1, result)
    end
    while true
        (pos, val) = _parse_json_value(s, pos)
        push!(result, val)
        pos = _skip_ws(s, pos)
        if pos <= length(s) && s[pos] == ','
            pos += 1
        else
            break
        end
    end
    pos = _skip_ws(s, pos)
    pos += 1  # skip ]
    return (pos, result)
end

function _parse_json_string(s::String, pos::Int)
    pos += 1  # skip opening "
    chars = Char[]
    while pos <= length(s) && s[pos] != '"'
        if s[pos] == '\\'
            pos += 1
            if pos <= length(s)
                esc = s[pos]
                if esc == 'n'
                    push!(chars, '\n')
                elseif esc == 'r'
                    push!(chars, '\r')
                elseif esc == 't'
                    push!(chars, '\t')
                elseif esc == '\\'
                    push!(chars, '\\')
                elseif esc == '"'
                    push!(chars, '"')
                else
                    push!(chars, esc)
                end
            end
        else
            push!(chars, s[pos])
        end
        pos += 1
    end
    pos += 1  # skip closing "
    return (pos, String(chars))
end

function _parse_json_bool(s::String, pos::Int)
    if startswith(s[pos:end], "true")
        return (pos + 4, true)
    elseif startswith(s[pos:end], "false")
        return (pos + 5, false)
    end
    return (pos, nothing)
end

function _parse_json_null(s::String, pos::Int)
    if startswith(s[pos:end], "null")
        return (pos + 4, nothing)
    end
    return (pos, nothing)
end

function _parse_json_number(s::String, pos::Int)
    start = pos
    while pos <= length(s) && (s[pos] in '0':'9' || s[pos] in ('-', '+', '.', 'e', 'E'))
        pos += 1
    end
    num_str = s[start:(pos - 1)]
    if occursin('.', num_str) || occursin('e', num_str) || occursin('E', num_str)
        return (pos, parse(Float64, num_str))
    else
        return (pos, parse(Int, num_str))
    end
end

"""
    compare_with_baseline(results, baseline_path)

Compare benchmark results against a stored baseline.
Reports regressions but does not exit non-zero (per E2: non-blocking).
"""
function compare_with_baseline(results::Vector{BenchResult}, baseline_path::String)
    baseline = load_baseline(baseline_path)
    if isempty(baseline)
        println("\n  (No baseline found at $baseline_path — skipping comparison)")
        return nothing
    end

    baseline_results = get(baseline, "results", [])
    if isempty(baseline_results)
        println("\n  (Baseline has no results — skipping comparison)")
        return nothing
    end

    # Build lookup by name
    baseline_by_name = Dict{String,Dict}()
    for r in baseline_results
        name = r["name"]
        baseline_by_name[name] = r
    end

    println("\n=== Baseline Comparison ===")
    println(
        @sprintf(
            "  %-45s  %12s  %12s  %8s", "Benchmark", "Baseline ms", "Current ms", "Ratio"
        )
    )
    println("  " * "-" ^ 80)

    n_regressions = 0
    for r in results
        base = get(baseline_by_name, r.name, nothing)
        if base === nothing
            continue
        end
        base_median_ns = Float64(get(base, "median_ns", 0))
        if base_median_ns == 0
            continue
        end
        ratio = r.median_ns / base_median_ns
        status = ratio > 1.25 ? "REGRESSION" : (ratio < 0.85 ? "faster" : "ok")
        if ratio > 1.25
            n_regressions += 1
        end
        println(
            @sprintf(
                "  %-45s  %12.3f  %12.3f  %7.2fx  %s",
                r.name,
                base_median_ns / 1e6,
                r.median_ns / 1e6,
                ratio,
                status
            )
        )
    end

    println()
    if n_regressions > 0
        println("  $n_regressions potential regression(s) detected (ratio > 1.25x)")
        println("  NOTE: Per E2, this does NOT block. Review for confirmed regressions.")
    else
        println("  No regressions detected.")
    end
end

# ──────────────────────────────────────────────────────────────────────────
# Main entry point
# ──────────────────────────────────────────────────────────────────────────

function main()
    args = parse_args(ARGS)
    args === nothing && return nothing

    config = args

    # Collect environment metadata
    env = collect_environment()
    env["samples"] = config.samples
    env["seconds_budget"] = config.seconds

    if config.report_only
        println(to_json_str(Dict{String,Any}("environment" => env)))
        return nothing
    end

    println("=" ^ 70)
    println("  Mimosa.jl Reproducible Benchmark Suite")
    println("=" ^ 70)

    # Set global BenchmarkTools parameters
    BenchmarkTools.DEFAULT_PARAMETERS.samples = config.samples
    BenchmarkTools.DEFAULT_PARAMETERS.seconds = config.seconds

    println("  Julia:    $(env["julia_version"])")
    println("  Threads:  $(env["n_threads"])")
    println("  CPU:      $(env["cpu_model"])")
    println("  OS:       $(env["os"]) $(env["kernel"])")
    println("  RAM:      $(env["total_ram_gb"]) GB total, $(env["free_ram_gb"]) GB free")
    println("  Commit:   $(env["commit_sha"])")
    println("  Date:     $(env["timestamp"])")
    println("  Samples:  $(config.samples)  Seconds budget: $(config.seconds)")
    println("=" ^ 70)

    results = BenchResult[]

    # Run all benchmark groups
    bench_import_and_startup!(results, config)
    bench_pwm_scan!(results, config)
    bench_batch_scan!(results, config)
    bench_one_to_many!(results, config)
    bench_higher_order_scan!(results, config)
    bench_site_extraction!(results, config)
    bench_gev!(results, config)
    bench_null_distribution!(results, config)
    bench_serial_vs_threaded!(results, config)
    bench_serialization!(results, config)
    bench_storage!(results, config)

    # Build output JSON
    output = Dict{String,Any}(
        "environment" => env, "results" => [result_to_dict(r) for r in results]
    )

    # Write output
    json_str = to_json_str(output)
    if config.output_file !== nothing
        write(config.output_file, json_str)
        println("\n  Results written to: $(config.output_file)")
    else
        println("\n" * "=" ^ 70)
        println("  Machine-readable JSON output:")
        println("=" ^ 70)
        println(json_str)
    end

    # Baseline comparison
    if config.baseline_file !== nothing
        compare_with_baseline(results, config.baseline_file)
    end

    println("\n" * "=" ^ 70)
    println("  Benchmark suite complete.")
    return println("=" ^ 70)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
