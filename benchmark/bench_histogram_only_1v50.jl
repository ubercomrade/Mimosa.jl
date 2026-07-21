#!/usr/bin/env julia
# Experimental 1-vs-50 benchmark for histogram-only normalization.
#
# Unlike HybridEmpiricalLogTail, HistogramOnlyEmpiricalLogTail never collects
# or sorts an exact high-score tail. The benchmark sweeps histogram sizes and
# compares speed and output differences against exact empirical normalization.
# This strategy is intentionally benchmark-local until its accuracy trade-off
# is understood.

include(joinpath(@__DIR__, "bench_1v50.jl"))

const N_BACKGROUND = parse(Int, get(ENV, "MIMOSA_BENCH_BG_SEQUENCES", "500"))
const BACKGROUND_LENGTH = parse(Int, get(ENV, "MIMOSA_BENCH_BG_LENGTH", "20000"))

struct HistogramOnlyEmpiricalLogTail <: AbstractNormalizationStrategy
    bins::Int

    function HistogramOnlyEmpiricalLogTail(bins::Integer)
        value = Int(bins)
        256 <= value <= 1_048_576 ||
            throw(ArgumentError("histogram bins must be in 256:1_048_576"))
        ispow2(value) || throw(ArgumentError("histogram bins must be a power of two"))
        return new(value)
    end
end

struct HistogramOnlyLogTailTable
    minimum::Float32
    bin_width::Float32
    log_tail::Vector{Float32}
end

function Mimosa.normalization_fingerprint(strategy::HistogramOnlyEmpiricalLogTail)
    return "histogram-only-log-tail-v1;bins=$(strategy.bins)"
end

@inline function histogram_index(
    value::Float32, minimum::Float32, bin_width::Float32, bins::Int
)
    bin_width == 0.0f0 && return 1
    return clamp(floor(Int, (value - minimum) / bin_width) + 1, 1, bins)
end

function fit_histogram_only(
    strategy::HistogramOnlyEmpiricalLogTail,
    scores::AbstractVector{T};
    scan_execution::ExecutionPolicy=SerialExecution(),
) where {T<:Real}
    values = Float32.(scores)
    all(isfinite, values) || throw(ArgumentError("normalization scores must be finite"))
    n = length(values)
    n == 0 && return HistogramOnlyLogTailTable(0.0f0, 1.0f0, Float32[0.0f0])

    minimum_score, maximum_score = extrema(values)
    bins = minimum_score == maximum_score ? 1 : strategy.bins
    bin_width = minimum_score == maximum_score ?
                0.0f0 : (maximum_score - minimum_score) / Float32(bins)
    nchunks = scan_execution isa ThreadedExecution ?
              Mimosa._effective_ntasks(scan_execution, n) : 1
    counts = zeros(UInt32, bins, nchunks)

    Mimosa._parallel_for(scan_execution, nchunks) do chunk
        first_index = fld((chunk - 1) * n, nchunks) + 1
        last_index = fld(chunk * n, nchunks)
        @inbounds for index in first_index:last_index
            bin = histogram_index(values[index], minimum_score, bin_width, bins)
            counts[bin, chunk] += UInt32(1)
        end
    end

    log_tail = Vector{Float32}(undef, bins)
    cumulative = UInt64(0)
    @inbounds for bin in bins:-1:1
        for chunk in 1:nchunks
            cumulative += counts[bin, chunk]
        end
        log_tail[bin] = Float32(-log10(Float64(cumulative) / Float64(n)))
    end
    return HistogramOnlyLogTailTable(minimum_score, bin_width, log_tail)
end

function transform_histogram_only(
    table::HistogramOnlyLogTailTable,
    scores::RaggedArray{Float32};
    scan_execution::ExecutionPolicy=SerialExecution(),
)
    n = length(scores.data)
    output = Vector{Float32}(undef, n)
    nchunks = scan_execution isa ThreadedExecution ?
              Mimosa._effective_ntasks(scan_execution, max(n, 1)) : 1
    Mimosa._parallel_for(scan_execution, nchunks) do chunk
        first_index = fld((chunk - 1) * n, nchunks) + 1
        last_index = fld(chunk * n, nchunks)
        @inbounds for index in first_index:last_index
            bin = histogram_index(
                scores.data[index], table.minimum, table.bin_width, length(table.log_tail)
            )
            output[index] = table.log_tail[bin]
        end
    end
    return RaggedArray(output, copy(scores.offsets))
end

function normalize_histogram_only(
    table::HistogramOnlyLogTailTable,
    bundle::StrandPair{<:RaggedArray{Float32}};
    scan_execution::ExecutionPolicy=SerialExecution(),
)
    forward = transform_histogram_only(
        table, bundle.forward; scan_execution=scan_execution
    )
    bundle.forward === bundle.reverse && return StrandPair(forward, forward)
    reverse = transform_histogram_only(
        table, bundle.reverse; scan_execution=scan_execution
    )
    return StrandPair(forward, reverse)
end

function Mimosa._fit_normalize(
    strategy::HistogramOnlyEmpiricalLogTail,
    raw::StrandPair{<:RaggedArray{Float32}};
    calibration::StrandPair{<:RaggedArray{Float32}}=raw,
    tail_logfpr::Real=0.0,
    scan_execution::ExecutionPolicy=SerialExecution(),
)
    isfinite(tail_logfpr) && tail_logfpr >= 0 ||
        throw(ArgumentError("tail_logfpr must be finite and non-negative"))
    table = fit_histogram_only(
        strategy,
        Mimosa._empirical_workspace(calibration);
        scan_execution=scan_execution,
    )
    normalized = normalize_histogram_only(
        table, raw; scan_execution=scan_execution
    )
    return table, normalized
end

function detailed_differences(exact_results, approximate_results)
    base = compare_outputs(exact_results, approximate_results)
    signed = getproperty.(approximate_results, :score) .- getproperty.(exact_results, :score)
    absolute = sort(abs.(signed))
    site_deltas = abs.(
        getproperty.(approximate_results, :n_sites) .- getproperty.(exact_results, :n_sites)
    )
    worst_index = argmax(abs.(signed))
    p95_index = clamp(ceil(Int, 0.95 * length(absolute)), 1, length(absolute))
    return merge(
        base,
        (
            mean_signed_score_delta=sum(signed) / length(signed),
            p95_score_delta=absolute[p95_index],
            maximum_site_delta=maximum(site_deltas),
            total_site_delta=sum(site_deltas),
            worst_target=exact_results[worst_index].target,
        ),
    )
end

function parse_histogram_bins()
    text = get(ENV, "MIMOSA_BENCH_HISTOGRAM_BINS", "65536,262144,1048576")
    values = parse.(Int, strip.(split(text, ',')))
    isempty(values) && error("MIMOSA_BENCH_HISTOGRAM_BINS must not be empty")
    for bins in values
        HistogramOnlyEmpiricalLogTail(bins)
    end
    return values
end

function show_summary_row(label::String, result, exact_results, exact_ns::Float64)
    differences = detailed_differences(exact_results, result.results)
    milliseconds = Float64(result.many.median_ns) / 1.0e6
    speedup = exact_ns / Float64(result.many.median_ns)
    println(
        rpad(label, 31),
        @sprintf(
            "%9.3f ms  %6.2fx  mean|Δ|=%-10.7f p95=%-10.7f max=%-10.7f  offsets=%d orientations=%d sites=%d maxΔsites=%d",
            milliseconds,
            speedup,
            differences.mean_score_delta,
            differences.p95_score_delta,
            differences.max_score_delta,
            differences.offset_changes,
            differences.orientation_changes,
            differences.site_changes,
            differences.maximum_site_delta,
        ),
    )
    println(
        "  " ^ 31,
        @sprintf(
            "score bias=%+.7f  worst target=%s  total |Δsites|=%d",
            differences.mean_signed_score_delta,
            differences.worst_target,
            differences.total_site_delta,
        ),
    )
    return differences
end

function main_histogram_only()
    histogram_bins = parse_histogram_bins()
    threaded = ThreadedExecution()
    exact_case = BenchmarkCase(
        "Exact target-threaded",
        Mimosa.EmpiricalLogTail(),
        threaded,
        SerialExecution(),
    )
    hybrid_case = BenchmarkCase(
        "Hybrid exact-tail",
        HybridEmpiricalLogTail(HYBRID_BINS),
        SerialExecution(),
        threaded,
    )
    histogram_cases = [
        BenchmarkCase(
            "Histogram-only $bins bins",
            HistogramOnlyEmpiricalLogTail(bins),
            SerialExecution(),
            threaded,
        ) for bins in histogram_bins
    ]

    println("=" ^ 108)
    println("  Mimosa.jl — histogram-only accuracy/speed benchmark")
    println("=" ^ 108)
    println("  Workload           : 1-vs-$N_TARGETS on $N_SEQUENCES × $SEQ_LENGTH bp")
    println("  Calibration        : $N_BACKGROUND × $BACKGROUND_LENGTH bp")
    println("  Minimum log FPR    : $MIN_LOGFPR")
    println("  Tail Hybrid bins   : $HYBRID_BINS")
    println("  Histogram-only bins: $(join(histogram_bins, ", "))")
    println("  Repetitions        : $N_REPS")
    println("  Julia threads      : $(Threads.nthreads())")
    println("  Date               : $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
    println("=" ^ 108)

    sequences = make_random_sequences(N_SEQUENCES, SEQ_LENGTH; seed=SEED)
    background = make_random_sequences(
        N_BACKGROUND, BACKGROUND_LENGTH; seed=SEED + 1
    )
    query_model = make_pwm(PWM_WIDTH, SEED)
    target_models = [make_pwm(PWM_WIDTH, SEED + 100 + index) for index in 1:N_TARGETS]

    exact = benchmark_case(
        exact_case,
        query_model,
        target_models,
        sequences,
        N_REPS;
        background=background,
    )
    hybrid = benchmark_case(
        hybrid_case,
        query_model,
        target_models,
        sequences,
        N_REPS;
        background=background,
    )
    histogram_results = [
        benchmark_case(
            case,
            query_model,
            target_models,
            sequences,
            N_REPS;
            background=background,
        ) for
        case in histogram_cases
    ]

    exact_ns = Float64(exact.many.median_ns)
    println("\n" * "=" ^ 108)
    println("  SUMMARY — speedup and computation differences relative to Exact")
    println("=" ^ 108)
    println(rpad("Exact target-threaded", 31), @sprintf("%9.3f ms  %6.2fx", exact_ns / 1.0e6, 1.0))
    show_summary_row("Hybrid exact-tail", hybrid, exact.results, exact_ns)
    for (case, result) in zip(histogram_cases, histogram_results)
        show_summary_row(case.label, result, exact.results, exact_ns)
    end
    println("=" ^ 108)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_histogram_only()
end
