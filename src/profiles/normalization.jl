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
struct EmpiricalLogTail end

function _fit_empirical_table!(workspace::Vector{Float32})
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
        value = Float32(-log10(Float64(j - 1) / Float64(n)))
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
    return _fit_empirical_table!(Float32.(scores))
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
    _fit_normalize_empirical(raw; calibration=raw, scan_execution=SerialExecution())

Canonical empirical-normalization pipeline. It fits a table only from
`calibration` and applies that table to `raw`.
"""
function _fit_normalize_empirical(
    raw::StrandPair{<:RaggedArray{Float32}};
    calibration::StrandPair{<:RaggedArray{Float32}}=raw,
    scan_execution::ExecutionPolicy=SerialExecution(),
)
    table = _fit_empirical_table!(_empirical_workspace(calibration))
    return table, normalize_bundle(table, raw; scan_execution=scan_execution)
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

"""
    transform_scores(table::LogTailTable, scores::RaggedArray{Float32})

Apply the log-tail lookup to every element of a [`RaggedArray`](@ref),
returning a new `RaggedArray{Float32}` with mapped values.
"""
function transform_scores(
    table::LogTailTable,
    scores::RaggedArray{Float32};
    scan_execution::ExecutionPolicy=SerialExecution(),
)
    n = length(scores.data)
    n == 0 && return RaggedArray(Float32[], copy(scores.offsets))
    out_data = Vector{Float32}(undef, n)
    nchunks = scan_execution isa ThreadedExecution ? _effective_ntasks(scan_execution, n) : 1
    _parallel_for(scan_execution, nchunks) do chunk
        first = fld((chunk - 1) * n, nchunks) + 1
        last = fld(chunk * n, nchunks)
        @inbounds for i in first:last
            out_data[i] = lookup_score(table, scores.data[i])
        end
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
    scan_execution::ExecutionPolicy=SerialExecution(),
)
    fwd = transform_scores(table, bundle.forward; scan_execution=scan_execution)
    bundle.forward === bundle.reverse && return StrandPair(fwd, fwd)
    rev = transform_scores(table, bundle.reverse; scan_execution=scan_execution)
    return StrandPair(fwd, rev)
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
