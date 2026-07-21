#!/usr/bin/env julia
# Benchmark normalization designs on the production-sized profile workload.
#
# It evaluates prototypes without changing Mimosa's public normalization API:
#   1. scan background directly into a contiguous calibration buffer;
#   2. exact lookup through a Float32 radix directory (no foreground sortperm);
#   3. one-pass parallel bucket partition followed by independent sorts;
#   5. reusable calibration workspace across targets;
#   7. approximate parallel histograms and their error against the exact table.
#
# Example:
#   JULIA_NUM_THREADS=4 julia --project=benchmark benchmark/bench_normalization_options.jl

using Printf
using Random
using Statistics
using Mimosa

const N_BACKGROUND = parse(Int, get(ENV, "MIMOSA_BENCH_BG_SEQUENCES", "500"))
const BACKGROUND_LENGTH = parse(Int, get(ENV, "MIMOSA_BENCH_BG_LENGTH", "20000"))
const N_FOREGROUND = parse(Int, get(ENV, "MIMOSA_BENCH_FG_SEQUENCES", "100"))
const FOREGROUND_LENGTH = parse(Int, get(ENV, "MIMOSA_BENCH_FG_LENGTH", "10000"))
const PWM_WIDTH = parse(Int, get(ENV, "MIMOSA_BENCH_PWM_WIDTH", "15"))
const N_REPEATS = parse(Int, get(ENV, "MIMOSA_BENCH_N_REPS", "3"))
const SEED = parse(Int, get(ENV, "MIMOSA_BENCH_SEED", "12345"))
const DIRECTORY_BITS = parse(Int, get(ENV, "MIMOSA_BENCH_DIRECTORY_BITS", "16"))
const PARTITION_BITS = parse(Int, get(ENV, "MIMOSA_BENCH_PARTITION_BITS", "12"))

function make_pwm(width::Int, seed::Int)
    rng = MersenneTwister(seed)
    weights = Matrix{Float32}(undef, 5, width)
    @inbounds for col in 1:width
        for row in 1:4
            weights[row, col] = Float32(randn(rng) * 0.5)
        end
        weights[5, col] = minimum(@view(weights[1:4, col])) - 0.1f0
    end
    return PWM("normalization_bench", weights, (0.25f0, 0.25f0, 0.25f0, 0.25f0))
end

function measure(f)
    f() # compile and allocate one-time method caches before measuring
    GC.gc()
    timed = @timed f()
    return (seconds=timed.time, bytes=timed.bytes, value=timed.value)
end

function median_measure(f)
    samples = [measure(f) for _ in 1:N_REPEATS]
    return samples[sortperm(getfield.(samples, :seconds))[cld(N_REPEATS, 2)]]
end

function median_timing(f)
    f() # compile outside the timed samples
    samples = NamedTuple{(:seconds, :bytes),Tuple{Float64,Int}}[]
    for _ in 1:N_REPEATS
        GC.gc()
        timed = @timed f()
        push!(samples, (seconds=timed.time, bytes=timed.bytes))
    end
    return samples[sortperm(getfield.(samples, :seconds))[cld(N_REPEATS, 2)]]
end

@inline function float_key(x::Float32)::UInt32
    bits = reinterpret(UInt32, x)
    return bits & 0x80000000 == 0 ? xor(bits, 0x80000000) : ~bits
end

struct RadixDirectory
    offsets::Vector{Int}
    bits::Int
end

function build_directory(scores::Vector{Float32}, bits::Int)
    1 <= bits <= 20 || throw(ArgumentError("directory bits must be in 1:20"))
    nbins = 1 << bits
    offsets = Vector{Int}(undef, nbins)
    index = 1
    shift = 32 - bits
    @inbounds for bin in (nbins - 1):-1:0
        while index <= length(scores) && Int(float_key(scores[index]) >> shift) > bin
            index += 1
        end
        offsets[bin + 1] = index
    end
    return RadixDirectory(offsets, bits)
end

@inline function lookup_directory(
    table::LogTailTable, directory::RadixDirectory, score::Float32
)::Float32
    scores = table.scores
    bin = Int(float_key(score) >> (32 - directory.bits))
    first = @inbounds directory.offsets[bin + 1]
    last = bin == 0 ? length(scores) : @inbounds(directory.offsets[bin]) - 1
    key = float_key(score)
    @inbounds while first <= last
        mid = (first + last) >>> 1
        if float_key(scores[mid]) > key
            first = mid + 1
        else
            last = mid - 1
        end
    end
    return table.log_tail[min(first, length(scores))]
end

function transform_directory!(
    output::Vector{Float32}, table::LogTailTable, directory::RadixDirectory, input::Vector{Float32}
)
    Threads.@threads :static for chunk in 1:Threads.nthreads()
        first = fld((chunk - 1) * length(input), Threads.nthreads()) + 1
        last = fld(chunk * length(input), Threads.nthreads())
        @inbounds for index in first:last
            output[index] = lookup_directory(table, directory, input[index])
        end
    end
    return output
end

function scan_calibration_direct(model::AbstractMotifModel, batch::EncodedSequenceBatch)
    offsets = Mimosa._scan_offsets(batch, model)
    n = offsets[end] - 1
    data = Vector{Float32}(undef, 2 * n)
    Threads.@threads :static for row_index in 1:nsequences(batch)
        first = offsets[row_index]
        last = offsets[row_index + 1] - 1
        length = last - first + 1
        forward = @view(data[first:last])
        reverse = @view(data[n + first:n + last])
        Mimosa.scan_both!(forward, reverse, model, sequence(batch, row_index), length)
    end
    return data
end

function scan_calibration_current(model::AbstractMotifModel, batch::EncodedSequenceBatch)
    raw = scan(model, batch; strands=BothStrands(), execution=ThreadedExecution())
    return Mimosa._empirical_workspace(raw)
end

function partition_sort!(destination::Vector{Float32}, source::Vector{Float32}, bits::Int)
    nworkers = Threads.nthreads()
    nbins = 1 << bits
    shift = 32 - bits
    counts = zeros(Int, nbins, nworkers)
    Threads.@threads :static for worker in 1:nworkers
        first = fld((worker - 1) * length(source), nworkers) + 1
        last = fld(worker * length(source), nworkers)
        @inbounds for index in first:last
            bin = Int(float_key(source[index]) >> shift) + 1
            counts[bin, worker] += 1
        end
    end

    positions = similar(counts)
    next = 1
    @inbounds for bin in nbins:-1:1
        for worker in 1:nworkers
            positions[bin, worker] = next
            next += counts[bin, worker]
        end
    end

    Threads.@threads :static for worker in 1:nworkers
        first = fld((worker - 1) * length(source), nworkers) + 1
        last = fld(worker * length(source), nworkers)
        local_positions = @view positions[:, worker]
        @inbounds for index in first:last
            value = source[index]
            bin = Int(float_key(value) >> shift) + 1
            destination[local_positions[bin]] = value
            local_positions[bin] += 1
        end
    end

    starts = Vector{Int}(undef, nbins)
    stops = Vector{Int}(undef, nbins)
    next = 1
    @inbounds for bin in nbins:-1:1
        starts[bin] = next
        next += sum(@view counts[bin, :])
        stops[bin] = next - 1
    end
    Threads.@threads for bin in 1:nbins
        start = starts[bin]
        stop = stops[bin]
        start < stop && sort!(@view(destination[start:stop]); rev=true)
    end
    return destination
end

mutable struct CalibrationWorkspace
    scores::Vector{Float32}
end

CalibrationWorkspace(n::Int) = CalibrationWorkspace(Vector{Float32}(undef, n))

function fit_reusing_workspace!(workspace::CalibrationWorkspace, source::Vector{Float32})
    resize!(workspace.scores, length(source))
    copyto!(workspace.scores, source)
    return Mimosa._fit_empirical_table!(workspace.scores)
end

struct HistogramTable
    lower::Float32
    width::Float32
    tail::Vector{Float32}
end

struct HybridHistogramTable
    histogram::HistogramTable
    exact_tail::LogTailTable
end

@inline function histogram_bin(score::Float32, lower::Float32, width::Float32, nbins::Int)
    return clamp(Int(floor((score - lower) / width)) + 1, 1, nbins)
end

function fit_histogram(scores::Vector{Float32}, nbins::Int)
    lower = Base.minimum(scores)
    upper = Base.maximum(scores)
    width = (upper - lower) / Float32(nbins)
    width > 0.0f0 || return HistogramTable(lower, 1.0f0, Float32[0.0f0])
    counts = zeros(Int, nbins, Threads.nthreads())
    Threads.@threads :static for worker in 1:Threads.nthreads()
        first = fld((worker - 1) * length(scores), Threads.nthreads()) + 1
        last = fld(worker * length(scores), Threads.nthreads())
        @inbounds for index in first:last
            bin = histogram_bin(scores[index], lower, width, nbins)
            counts[bin, worker] += 1
        end
    end
    tail = Vector{Float32}(undef, nbins)
    cumulative = 0
    @inbounds for bin in nbins:-1:1
        cumulative += sum(@view counts[bin, :])
        tail[bin] = Float32(-log10(Float64(cumulative) / Float64(length(scores))))
    end
    return HistogramTable(lower, width, tail)
end

@inline function lookup_histogram(table::HistogramTable, score::Float32)::Float32
    return table.tail[histogram_bin(score, table.lower, table.width, length(table.tail))]
end

@inline function lookup_hybrid(table::HybridHistogramTable, score::Float32)::Float32
    score >= table.exact_tail.scores[end] && return lookup_score(table.exact_tail, score)
    return lookup_histogram(table.histogram, score)
end

function approximation_error(
    exact::LogTailTable, approximate::HistogramTable, values::Vector{Float32}
)
    errors = Vector{Float32}(undef, length(values))
    tail_errors = Float32[]
    @inbounds for index in eachindex(values)
        reference = lookup_score(exact, values[index])
        error = abs(reference - lookup_histogram(approximate, values[index]))
        errors[index] = error
        reference >= 3.0f0 && push!(tail_errors, error)
    end
    sort!(errors)
    sort!(tail_errors)
    return (
        mean=mean(errors),
        p95=errors[cld(95 * length(errors), 100)],
        maximum=last(errors),
        tail_count=length(tail_errors),
        tail_p95=isempty(tail_errors) ? 0.0f0 : tail_errors[cld(95 * length(tail_errors), 100)],
        tail_max=isempty(tail_errors) ? 0.0f0 : last(tail_errors),
    )
end

function hybrid_error(
    exact::LogTailTable, approximate::HistogramTable, values::Vector{Float32}; tail_logfpr::Float32=3.0f0
)
    cutoff_index = findlast(>=(tail_logfpr), exact.log_tail)
    cutoff_index === nothing && return nothing
    hybrid = HybridHistogramTable(
        approximate,
        LogTailTable(exact.scores[1:cutoff_index], exact.log_tail[1:cutoff_index]),
    )
    errors = Vector{Float32}(undef, length(values))
    @inbounds for index in eachindex(values)
        errors[index] = abs(lookup_score(exact, values[index]) - lookup_hybrid(hybrid, values[index]))
    end
    sort!(errors)
    return (
        cutoff=hybrid.exact_tail.scores[end],
        entries=length(hybrid.exact_tail.scores),
        mean=mean(errors),
        p95=errors[cld(95 * length(errors), 100)],
        maximum=last(errors),
    )
end

function fit_sorted_table!(scores::Vector{Float32}, total::Int=length(scores))
    n = length(scores)
    n == 0 && return LogTailTable(Float32[0.0f0], Float32[0.0f0])
    log_tail = Vector{Float32}(undef, n)
    n_unique = 0
    k = 1
    @inbounds while k <= n
        score = scores[k]
        n_unique += 1
        scores[n_unique] = score
        j = k + 1
        while j <= n && scores[j] == score
            j += 1
        end
        log_tail[n_unique] = Float32(-log10(Float64(j - 1) / Float64(total)))
        k = j
    end
    resize!(scores, n_unique)
    resize!(log_tail, n_unique)
    return LogTailTable(scores, log_tail)
end

function collect_histogram_tail(scores::Vector{Float32}, histogram::HistogramTable, tail_logfpr::Float32)
    bin = findfirst(>=(tail_logfpr), histogram.tail)
    bin === nothing && return Float32[]
    cutoff = histogram.lower + Float32(bin - 1) * histogram.width
    tail = Float32[]
    sizehint!(tail, cld(length(scores), 1_000))
    @inbounds for score in scores
        score >= cutoff && push!(tail, score)
    end
    sort!(tail; rev=true)
    return tail
end

function normalize_hybrid_bundle(
    table::HybridHistogramTable,
    bundle::StrandPair{<:RaggedArray{Float32}},
)
    function transform(data::Vector{Float32})
        output = similar(data)
        Threads.@threads :static for worker in 1:Threads.nthreads()
            first = fld((worker - 1) * length(data), Threads.nthreads()) + 1
            last = fld(worker * length(data), Threads.nthreads())
            @inbounds for index in first:last
                output[index] = lookup_hybrid(table, data[index])
            end
        end
        return output
    end
    forward = RaggedArray(transform(bundle.forward.data), copy(bundle.forward.offsets))
    reverse = RaggedArray(transform(bundle.reverse.data), copy(bundle.reverse.offsets))
    return StrandPair(forward, reverse)
end

function prepare_exact_variant(
    model::AbstractMotifModel,
    foreground::StrandPair{<:RaggedArray{Float32}},
    background::EncodedSequenceBatch,
)
    calibration = scan_calibration_direct(model, background)
    sorted = similar(calibration)
    partition_sort!(sorted, calibration, PARTITION_BITS)
    table = fit_sorted_table!(sorted)
    return (table=table, bundle=normalize_bundle(table, foreground; scan_execution=ThreadedExecution()))
end

function prepare_hybrid_variant(
    model::AbstractMotifModel,
    foreground::StrandPair{<:RaggedArray{Float32}},
    background::EncodedSequenceBatch;
    bins::Int,
    tail_logfpr::Float32,
)
    calibration = scan_calibration_direct(model, background)
    histogram = fit_histogram(calibration, bins)
    exact_tail = collect_histogram_tail(calibration, histogram, tail_logfpr)
    tail_table = fit_sorted_table!(exact_tail, length(calibration))
    table = HybridHistogramTable(histogram, tail_table)
    return (table=table, bundle=normalize_hybrid_bundle(table, foreground))
end

function anchor_overlap(first::AnchorCSR, second::AnchorCSR)
    matched = 0
    @inbounds for row in 1:(length(first.offsets) - 1)
        i = first.offsets[row]
        j = second.offsets[row]
        i_stop = first.offsets[row + 1] - 1
        j_stop = second.offsets[row + 1] - 1
        while i <= i_stop && j <= j_stop
            left = first.positions[i]
            right = second.positions[j]
            if left == right
                matched += 1
                i += 1
                j += 1
            elseif left < right
                i += 1
            else
                j += 1
            end
        end
    end
    return matched
end

function compare_variants(
    exact_bundle::StrandPair{<:RaggedArray{Float32}},
    hybrid_bundle::StrandPair{<:RaggedArray{Float32}},
    query::PreparedProfile,
    threshold::Float32,
)
    exact_anchors = Mimosa._collect_both_anchors(exact_bundle, threshold)
    hybrid_anchors = Mimosa._collect_both_anchors(hybrid_bundle, threshold)
    config = Mimosa.ProfileConfig(
        metric=Mimosa.OverlapCoefficient(), search_range=10, window_radius=5,
        min_logfpr=threshold,
    )
    exact_result = Mimosa.profile_compare(
        query.bundle, query.anchors, exact_bundle, exact_anchors, config,
    )
    hybrid_result = Mimosa.profile_compare(
        query.bundle, query.anchors, hybrid_bundle, hybrid_anchors, config,
    )
    return (
        exact_result=exact_result,
        hybrid_result=hybrid_result,
        exact_anchors=(length(exact_anchors[1].positions), length(exact_anchors[2].positions)),
        hybrid_anchors=(length(hybrid_anchors[1].positions), length(hybrid_anchors[2].positions)),
        overlap=(anchor_overlap(exact_anchors[1], hybrid_anchors[1]), anchor_overlap(exact_anchors[2], hybrid_anchors[2])),
    )
end

function show_measure(label::String, result)
    println(
        rpad(label, 54),
        @sprintf("%8.3f ms  %8.1f MiB allocated", result.seconds * 1e3, result.bytes / 2.0^20),
    )
end

function main()
    println("Mimosa normalization options benchmark")
    println("  Julia threads: $(Threads.nthreads())")
    println("  background: $N_BACKGROUND × $BACKGROUND_LENGTH; foreground: $N_FOREGROUND × $FOREGROUND_LENGTH")
    model = make_pwm(PWM_WIDTH, SEED)
    background = make_random_sequences(N_BACKGROUND, BACKGROUND_LENGTH; seed=SEED)
    foreground = make_random_sequences(N_FOREGROUND, FOREGROUND_LENGTH; seed=SEED + 1)

    current = scan_calibration_current(model, background)
    direct = scan_calibration_direct(model, background)
    @assert current == direct
    println("\n[1] Direct contiguous calibration scan")
    show_measure("current: two strand arrays + calibration copy", median_measure(() -> scan_calibration_current(model, background)))
    show_measure("direct: one calibration buffer", median_measure(() -> scan_calibration_direct(model, background)))

    exact_table = Mimosa._fit_empirical_table!(copy(direct))
    foreground_raw = scan(model, foreground; strands=BothStrands(), execution=ThreadedExecution())
    foreground_values = vcat(foreground_raw.forward.data, foreground_raw.reverse.data)
    directory = build_directory(exact_table.scores, DIRECTORY_BITS)
    directory_output = similar(foreground_values)
    @inbounds for score in foreground_values
        @assert lookup_score(exact_table, score) == lookup_directory(exact_table, directory, score)
    end
    println("\n[2] Exact foreground lookup")
    show_measure("current: sortperm + linear merge (both strands)", median_measure(() -> normalize_bundle(exact_table, foreground_raw; scan_execution=ThreadedExecution())))
    show_measure("directory: independent exact lookups", median_measure(() -> transform_directory!(directory_output, exact_table, directory, foreground_values)))

    sorted = copy(direct)
    partitioned = similar(direct)
    partition_sort!(partitioned, sorted, PARTITION_BITS)
    sort!(sorted; rev=true)
    @assert sorted == partitioned
    partition_table = fit_sorted_table!(copy(partitioned))
    @assert partition_table.scores == exact_table.scores
    @assert partition_table.log_tail == exact_table.log_tail
    println("\n[3] Exact calibration sort")
    show_measure("current: Base.sort!", median_measure(() -> sort!(copy(direct); rev=true)))
    show_measure("partition + bucket-local sort", median_measure(() -> partition_sort!(partitioned, direct, PARTITION_BITS)))

    workspace = CalibrationWorkspace(length(direct))
    fit_reusing_workspace!(workspace, direct)
    println("\n[5] Reusable calibration workspace")
    show_measure("new workspace per fit", median_measure(() -> Mimosa._fit_empirical_table!(copy(direct))))
    show_measure("reuse score workspace", median_measure(() -> fit_reusing_workspace!(workspace, direct)))

    println("\n[7] Approximate histogram table; foreground error vs exact")
    for nbins in (256, 4096, 65536)
        measured = median_measure(() -> fit_histogram(direct, nbins))
        histogram = measured.value
        error = approximation_error(exact_table, histogram, foreground_values)
        show_measure("histogram ($nbins bins)", measured)
        println(
            @sprintf(
                "  error: mean=%.5f p95=%.5f max=%.5f | tail>=3: n=%d p95=%.5f max=%.5f",
                error.mean, error.p95, error.maximum, error.tail_count, error.tail_p95, error.tail_max,
            ),
        )
        hybrid = hybrid_error(exact_table, histogram, foreground_values)
        println(
            @sprintf(
                "  hybrid exact tail>=3: %d table entries, cutoff=%.5f; mean=%.5f p95=%.5f max=%.5f",
                hybrid.entries, hybrid.cutoff, hybrid.mean, hybrid.p95, hybrid.maximum,
            ),
        )
    end

    # The diagnostic sections intentionally keep several full calibration tables
    # alive for equality/error checks. Release them before measuring the composed
    # pipelines, otherwise their peak memory obscures the production-like case.
    current = nothing
    direct = nothing
    exact_table = nothing
    foreground_values = nothing
    directory_output = nothing
    sorted = nothing
    partitioned = nothing
    partition_table = nothing
    workspace = nothing
    GC.gc()

    println("\n[full] Composed exact (1+3+5) versus hybrid (1+5+7)")
    target_raw = foreground_raw
    exact_timing = median_timing(() -> prepare_exact_variant(model, target_raw, background))
    show_measure("exact: direct scan + partition sort + reuse-ready table", exact_timing)
    exact_bundle = prepare_exact_variant(model, target_raw, background).bundle
    GC.gc()
    for bins in (4096, 65536)
        hybrid_timing = median_timing(
            () -> prepare_hybrid_variant(
                model, target_raw, background; bins=bins, tail_logfpr=3.0f0,
            ),
        )
        show_measure("hybrid: direct scan + $bins-bin histogram + exact tail", hybrid_timing)
        hybrid_prepared = prepare_hybrid_variant(
            model, target_raw, background; bins=bins, tail_logfpr=3.0f0,
        )
        query_model = make_pwm(PWM_WIDTH, SEED + 99)
        for threshold in (0.0f0, 3.0f0)
            query = prepare_profile(
                query_model, foreground; background=background, min_logfpr=threshold,
                scan_execution=ThreadedExecution(),
            )
            compared = compare_variants(exact_bundle, hybrid_prepared.bundle, query, threshold)
            exact_result = compared.exact_result
            hybrid_result = compared.hybrid_result
            println(
                @sprintf(
                    "  bins=%d threshold=%.0f: anchors exact=(%d,%d), hybrid=(%d,%d), shared=(%d,%d); ",
                    bins, threshold,
                    compared.exact_anchors..., compared.hybrid_anchors..., compared.overlap...,
                ) * @sprintf(
                    "compare score %.7f→%.7f (Δ=%.7f), offset %d→%d, orient %s→%s, sites %d→%d",
                    exact_result[1], hybrid_result[1], hybrid_result[1] - exact_result[1],
                    exact_result[2], hybrid_result[2], exact_result[3], hybrid_result[3], exact_result[4], hybrid_result[4],
                ),
            )
        end
    end
end

main()
