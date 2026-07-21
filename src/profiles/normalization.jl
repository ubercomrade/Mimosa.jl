# Empirical log-tail normalization for score profiles.

"""
    LogTailTable

Empirical score-to-log-tail lookup table.

Fields:
- `scores::Vector{Float32}`: unique scores in descending order.
- `log_tail::Vector{Float32}`: corresponding `-log10(tail_probability)` values.
"""
struct LogTailTable
    scores::Vector{Float32}
    log_tail::Vector{Float32}
end

"""
    EmpiricalLogTail

Empirical `-log10(tail)` normalization strategy.
Fit from a flat sample of scores, then transform score profiles through
the lookup table. Matches Python's `build_score_log_tail_table`.
"""
abstract type AbstractNormalizationStrategy end

struct EmpiricalLogTail <: AbstractNormalizationStrategy end

"""
    HybridEmpiricalLogTail(; bins=65_536)

Default empirical `-log10(tail)` normalization strategy. The bulk of the
score distribution is represented by an equal-width histogram, while the
tail used for threshold anchors is retained exactly.
"""
struct HybridEmpiricalLogTail <: AbstractNormalizationStrategy
    bins::Int

    function HybridEmpiricalLogTail(bins::Integer)
        b = Int(bins)
        b >= 256 && b <= 1_048_576 ||
            throw(ArgumentError("hybrid normalization bins must be in 256:1_048_576."))
        ispow2(b) ||
            throw(ArgumentError("hybrid normalization bins must be a power of two."))
        return new(b)
    end
end

HybridEmpiricalLogTail(; bins::Integer=65_536) = HybridEmpiricalLogTail(bins)

"""
    normalization_fingerprint(strategy)

Return the stable identifier used to record normalization compatibility in
caches and null distributions.
"""
function normalization_fingerprint end

function normalization_fingerprint(::EmpiricalLogTail)
    return "empirical-log-tail-v1"
end

function normalization_fingerprint(s::HybridEmpiricalLogTail)
    return "hybrid-log-tail-v2;bins=$(s.bins)"
end

function _fit_empirical_table!(workspace::Vector{Float32}; total_n::Int=length(workspace))
    n = length(workspace)
    n == 0 && return LogTailTable(Float32[0.0f0], Float32[0.0f0])

    sort!(workspace; rev=true)
    log_tail = Vector{Float32}(undef, n)
    n_unique = 0
    k = 1
    @inbounds while k <= n
        score = workspace[k]
        n_unique += 1
        workspace[n_unique] = score
        j = k + 1
        while j <= n && workspace[j] == score
            j += 1
        end
        value = Float32(-log10(Float64(j - 1) / Float64(total_n)))
        log_tail[n_unique] = value
        k = j
    end
    resize!(workspace, n_unique)
    resize!(log_tail, n_unique)
    return LogTailTable(workspace, log_tail)
end

"""
    fit(::EmpiricalLogTail, scores::AbstractVector)

Build a [`LogTailTable`](@ref) from a flat sample of scores.

Algorithm (matching Python's `build_score_log_tail_table`):
1. Sort scores descending.
2. Extract unique scores with counts (descending order).
3. Compute cumulative counts → `tail_probability = cumcount / total`.
4. Transform to `-log10(tail_probability)`.
"""
function fit(::EmpiricalLogTail, scores::AbstractVector{T}) where {T<:Real}
    values = Float32.(scores)
    all(isfinite, values) || throw(ArgumentError("normalization scores must be finite."))
    return _fit_empirical_table!(values)
end

"""
    HybridLogTailTable

Fitted lookup table produced by [`HybridEmpiricalLogTail`](@ref). It combines
histogram values for the body of the distribution with an exact empirical
table for the selected tail.
"""
struct HybridLogTailTable
    minimum::Float32
    bin_width::Float32
    log_tail::Vector{Float32}
    exact_tail::LogTailTable
end

mutable struct HybridNormalizationWorkspace
    counts::Matrix{UInt32}
    tail_counts::Vector{Int}
    tail_offsets::Vector{Int}
    exact_tail::Vector{Float32}
end

function HybridNormalizationWorkspace(bins::Int, nchunks::Int)
    return HybridNormalizationWorkspace(
        zeros(UInt32, bins, nchunks),
        zeros(Int, nchunks),
        zeros(Int, nchunks + 1),
        Float32[],
    )
end

function _copy_score_summary(scores::AbstractVector{<:Real}, execution::Execution)
    n = length(scores)
    values = Vector{Float32}(undef, n)
    n == 0 && return values, 0.0f0, 0.0f0
    nchunks = _effective_ntasks(execution, n)
    minima = fill(Inf32, nchunks)
    maxima = fill(-Inf32, nchunks)
    source_first = firstindex(scores)
    _parallel_chunks(execution, n) do first, last, chunk
        local_minimum = Inf32
        local_maximum = -Inf32
        finite = true
        @inbounds for i in first:last
            value = Float32(scores[source_first + i - 1])
            values[i] = value
            finite &= isfinite(value)
            local_minimum = min(local_minimum, value)
            local_maximum = max(local_maximum, value)
        end
        minima[chunk] = finite ? local_minimum : NaN32
        return maxima[chunk] = finite ? local_maximum : NaN32
    end
    all(isfinite, minima) && all(isfinite, maxima) ||
        throw(ArgumentError("normalization scores must be finite."))
    return values, minimum(minima), maximum(maxima)
end

function fit(
    strategy::HybridEmpiricalLogTail,
    scores::AbstractVector{T};
    tail_logfpr::Real=0.0,
    execution::Execution=Execution(),
) where {T<:Real}
    values, lo, hi = _copy_score_summary(scores, execution)
    n = length(scores)
    n == 0 && return HybridLogTailTable(
        0.0f0, 1.0f0, Float32[0.0f0], LogTailTable(Float32[0.0f0], Float32[0.0f0])
    )
    width = lo == hi ? 1.0f0 : (hi - lo) / Float32(strategy.bins)
    bins = lo == hi ? 1 : strategy.bins
    nchunks = _effective_ntasks(execution, n)
    ws = HybridNormalizationWorkspace(bins, nchunks)
    _parallel_chunks(execution, n) do first, last, chunk
        if first <= last
            @inbounds for i in first:last
                index = if lo == hi
                    1
                else
                    min(bins, max(1, floor(Int, (values[i] - lo) / width) + 1))
                end
                ws.counts[index, chunk] += UInt32(1)
            end
        end
    end
    counts = vec(sum(ws.counts; dims=2))
    cumulative = zeros(UInt64, bins)
    running = UInt64(0)
    @inbounds for i in bins:-1:1
        running += counts[i]
        cumulative[i] = running
    end
    isfinite(tail_logfpr) || throw(ArgumentError("tail_logfpr must be finite."))
    effective_tail_logfpr = max(0.0, Float64(tail_logfpr))
    cutoff_count = max(
        UInt64(1), ceil(UInt64, Float64(n) * 10.0 ^ (-effective_tail_logfpr))
    )
    # `cumulative` is indexed from low to high scores and decreases with the
    # index; choose the highest bin whose tail still reaches the threshold.
    cutoff_bin = findlast(>=(cutoff_count), cumulative)
    cutoff_bin = something(cutoff_bin, 1)

    _parallel_chunks(execution, n) do first, last, chunk
        count = 0
        @inbounds for i in first:last
            index =
                lo == hi ? 1 : min(bins, max(1, floor(Int, (values[i] - lo) / width) + 1))
            count += index >= cutoff_bin
        end
        return ws.tail_counts[chunk] = count
    end
    ws.tail_offsets[1] = 1
    @inbounds for chunk in 1:nchunks
        ws.tail_offsets[chunk + 1] = ws.tail_offsets[chunk] + ws.tail_counts[chunk]
    end
    resize!(ws.exact_tail, ws.tail_offsets[end] - 1)
    _parallel_chunks(execution, n) do first, last, chunk
        destination = ws.tail_offsets[chunk]
        @inbounds for i in first:last
            index =
                lo == hi ? 1 : min(bins, max(1, floor(Int, (values[i] - lo) / width) + 1))
            if index >= cutoff_bin
                ws.exact_tail[destination] = values[i]
                destination += 1
            end
        end
    end
    exact = _fit_empirical_table!(ws.exact_tail; total_n=n)
    histogram_log_tail = Vector{Float32}(undef, bins)
    @inbounds for i in 1:bins
        histogram_log_tail[i] = Float32(-log10(Float64(cumulative[i]) / n))
    end
    return HybridLogTailTable(lo, width, histogram_log_tail, exact)
end

function _empirical_workspace(bundle::StrandPair{<:RaggedArray{Float32}})
    fwd = bundle.forward.data
    rev = bundle.reverse.data
    n_fwd = length(fwd)
    n_rev = bundle.forward === bundle.reverse ? 0 : length(rev)
    total = Base.checked_add(n_fwd, n_rev)
    workspace = Vector{Float32}(undef, total)
    copyto!(workspace, 1, fwd, 1, n_fwd)
    n_rev > 0 && copyto!(workspace, n_fwd + 1, rev, 1, n_rev)
    return workspace
end

"""
    _fit_normalize_empirical(raw; calibration=raw, execution=Execution())

Canonical empirical-normalization pipeline. It fits a table only from
`calibration` and applies that table to `raw`.
"""
function _fit_normalize_empirical(
    raw::StrandPair{<:RaggedArray{Float32}};
    calibration::StrandPair{<:RaggedArray{Float32}}=raw,
    execution::Execution=Execution(),
)
    table = _fit_empirical_table!(_empirical_workspace(calibration))
    return table, normalize_bundle(table, raw; execution=execution)
end

function _fit_normalize(
    strategy::EmpiricalLogTail,
    raw::StrandPair{<:RaggedArray{Float32}};
    calibration::StrandPair{<:RaggedArray{Float32}}=raw,
    tail_logfpr::Real=0.0,
    execution::Execution=Execution(),
)
    return _fit_normalize_empirical(raw; calibration, execution)
end

function _fit_normalize(
    strategy::HybridEmpiricalLogTail,
    raw::StrandPair{<:RaggedArray{Float32}};
    calibration::StrandPair{<:RaggedArray{Float32}}=raw,
    tail_logfpr::Real=0.0,
    execution::Execution=Execution(),
)
    table = fit(strategy, _empirical_workspace(calibration); tail_logfpr, execution)
    return table, normalize_bundle(table, raw; execution)
end

"""
    _lower_bound_desc(scores::Vector{Float32}, target::Float32)

Binary search in a descending-sorted array: return the 1-based index of
the first element `<= target`. Matches Python's `_lower_bound_desc`.

Edge cases:
- `n <= 1`: return 1.
- `target >= scores[1]` (largest): return 1.
- `target <= scores[end]` (smallest): return n.
"""
function _lower_bound_desc(scores::Vector{Float32}, target::Float32)
    n = length(scores)
    n <= 1 && return 1
    @inbounds if target >= scores[1]
        return 1
    end
    @inbounds if target <= scores[n]
        return n
    end
    lo, hi = 1, n + 1
    @inbounds while lo < hi
        mid = (lo + hi) ÷ 2
        if scores[mid] > target
            lo = mid + 1
        else
            hi = mid
        end
    end
    return lo
end

"""
    lookup_score(table::LogTailTable, score::Float32)

Map a single score to its empirical log-tail value via descending binary search.
"""
function lookup_score(table::LogTailTable, score::Float32)
    idx = _lower_bound_desc(table.scores, score)
    return table.log_tail[idx]
end

function lookup_score(table::HybridLogTailTable, score::Float32)
    if !isempty(table.exact_tail.scores) && score >= table.exact_tail.scores[end]
        return lookup_score(table.exact_tail, score)
    end
    isempty(table.log_tail) && return 0.0f0
    index =
        table.bin_width == 0 ? 1 : floor(Int, (score - table.minimum) / table.bin_width) + 1
    return table.log_tail[clamp(index, 1, length(table.log_tail))]
end

function transform_scores(
    table::HybridLogTailTable,
    scores::RaggedArray{Float32};
    execution::Execution=Execution(),
)
    n = length(scores.data)
    out = Vector{Float32}(undef, n)
    _parallel_chunks(execution, n) do first, last, _
        @inbounds for i in first:last
            out[i] = lookup_score(table, scores.data[i])
        end
    end
    return RaggedArray(out, copy(scores.offsets))
end

"""
    _transform_scores_sorted!(out_data, table, input, first, last)

Transform one contiguous input range exactly as [`lookup_score`](@ref), while
accessing the large descending calibration table sequentially. Input indices
are sorted by score, then each score advances a single table cursor to the
first calibration score `<=` it. Results are written back in original order.
"""
function _transform_scores_sorted!(
    out_data::Vector{Float32},
    table::LogTailTable,
    input::Vector{Float32},
    first::Int,
    last::Int,
)
    n = last - first + 1
    n <= 0 && return nothing

    permutation = sortperm(@view(input[first:last]); rev=true)
    table_index = 1
    table_scores = table.scores
    table_log_tail = table.log_tail
    n_table = length(table_scores)

    @inbounds for local_index in permutation
        input_index = first + local_index - 1
        score = input[input_index]
        while table_index < n_table && table_scores[table_index] > score
            table_index += 1
        end
        out_data[input_index] = table_log_tail[table_index]
    end
    return nothing
end

"""
    transform_scores(table::LogTailTable, scores::RaggedArray{Float32})

Apply empirical log-tail normalization to every score. Each execution chunk
sorts its input scores descending and performs a linear merge against the
descending lookup table. This is exactly equivalent to applying
[`lookup_score`](@ref) individually, but avoids cache-unfriendly binary
searches in a large calibration table.
"""
function transform_scores(
    table::LogTailTable, scores::RaggedArray{Float32}; execution::Execution=Execution()
)
    n = length(scores.data)
    n == 0 && return RaggedArray(Float32[], copy(scores.offsets))
    out_data = Vector{Float32}(undef, n)
    _parallel_chunks(execution, n) do first, last, _
        return _transform_scores_sorted!(out_data, table, scores.data, first, last)
    end
    return RaggedArray(out_data, copy(scores.offsets))
end

"""
    flatten_bundle(bundle::StrandPair{<:RaggedArray{Float32}})

Flatten all valid scores from both strands into a single vector.
Used to fit the normalization table from background scan scores.
"""
function flatten_bundle(bundle::StrandPair{<:RaggedArray{Float32}})
    fwd = bundle.forward.data
    rev = bundle.reverse.data
    workspace = Vector{Float32}(undef, Base.checked_add(length(fwd), length(rev)))
    copyto!(workspace, 1, fwd, 1, length(fwd))
    copyto!(workspace, length(fwd) + 1, rev, 1, length(rev))
    return workspace
end

"""
    normalize_bundle(table::LogTailTable, bundle::StrandPair{<:RaggedArray{Float32}})

Apply the log-tail lookup to both strands of a profile bundle.
"""
function normalize_bundle(
    table::LogTailTable,
    bundle::StrandPair{<:RaggedArray{Float32}};
    execution::Execution=Execution(),
)
    fwd = transform_scores(table, bundle.forward; execution=execution)
    bundle.forward === bundle.reverse && return StrandPair(fwd, fwd)
    rev = transform_scores(table, bundle.reverse; execution=execution)
    return StrandPair(fwd, rev)
end

function normalize_bundle(
    table::HybridLogTailTable,
    bundle::StrandPair{<:RaggedArray{Float32}};
    execution::Execution=Execution(),
)
    fwd = transform_scores(table, bundle.forward; execution)
    bundle.forward === bundle.reverse && return StrandPair(fwd, fwd)
    return StrandPair(fwd, transform_scores(table, bundle.reverse; execution))
end

"""
    lookup_score_for_tail_probability(table::LogTailTable, tail_probability::Float64)

Convert a tail probability threshold to the corresponding score cutoff.
Matches Python's `lookup_score_for_tail_probability`.
"""
function lookup_score_for_tail_probability(table::LogTailTable, tail_probability::Float64)
    if tail_probability <= 0
        return Float64(table.scores[1])
    end
    target_log_tail = Float64(-log10(tail_probability))
    n = length(table.log_tail)
    # Find the last index where log_tail >= target (descending scores → ascending log_tail)
    last_valid = 0
    @inbounds for i in 1:n
        if Float64(table.log_tail[i]) >= target_log_tail
            last_valid = i
        end
    end
    if last_valid == 0
        return Float64(table.scores[end])
    end
    return Float64(table.scores[last_valid])
end
