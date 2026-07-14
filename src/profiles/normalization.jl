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
    n = length(scores)
    n == 0 && return LogTailTable(Float32[0.0f0], Float32[0.0f0])

    # Collect and sort descending
    flat = Float32[float(s) for s in scores]
    sort!(flat; rev=true)

    # Count unique scores in descending order.
    n_unique = 1
    @inbounds for i in 2:n
        if flat[i] != flat[i - 1]
            n_unique += 1
        end
    end

    unique_scores = Vector{Float32}(undef, n_unique)
    unique_scores[1] = flat[1]
    j = 1
    @inbounds for i in 2:n
        if flat[i] != unique_scores[j]
            j += 1
            unique_scores[j] = flat[i]
        end
    end

    # Cumulative counts → tail probabilities → -log10
    log_tail = Vector{Float32}(undef, n_unique)
    cum = 0
    group_start = 1
    @inbounds for i in 1:n_unique
        while group_start <= n && flat[group_start] == unique_scores[i]
            cum += 1
            group_start += 1
        end
        log_tail[i] = Float32(-log10(Float64(cum) / Float64(n)))
    end

    return LogTailTable(unique_scores, log_tail)
end

function _fit_transform_empirical(scores::RaggedArray{Float32})
    n = length(scores.data)
    n == 0 && return (
        LogTailTable(Float32[0.0f0], Float32[0.0f0]),
        RaggedArray(Float32[], copy(scores.offsets)),
    )

    perm = sortperm(scores.data; rev=true)
    n_unique = 1
    @inbounds for k in 2:n
        n_unique += scores.data[perm[k]] != scores.data[perm[k - 1]]
    end
    unique_scores = Vector{Float32}(undef, n_unique)
    log_tail = Vector{Float32}(undef, n_unique)
    normalized = Vector{Float32}(undef, n)
    group = 1
    cum = 0
    k = 1
    @inbounds while k <= n
        score = scores.data[perm[k]]
        unique_scores[group] = score
        j = k
        while j <= n && scores.data[perm[j]] == score
            j += 1
        end
        cum = j - 1
        value = Float32(-log10(Float64(cum) / Float64(n)))
        log_tail[group] = value
        for p in k:(j - 1)
            normalized[perm[p]] = value
        end
        group += 1
        k = j
    end
    table = LogTailTable(unique_scores, log_tail)
    return table, RaggedArray(normalized, copy(scores.offsets))
end

"""
    _fit_transform_empirical(bundle::StrandPair)

Fit one empirical table from both strands and scatter normalized values while
walking the single descending permutation. This is the common case where the
calibration bundle and the output bundle are the same raw scan.
"""
function _fit_transform_empirical(bundle::StrandPair{<:RaggedArray{Float32}})
    if bundle.forward === bundle.reverse
        table, normalized = _fit_transform_empirical(bundle.forward)
        return table, StrandPair(normalized, normalized)
    end

    n_forward = length(bundle.forward.data)
    n_reverse = length(bundle.reverse.data)
    n = n_forward + n_reverse
    if n == 0
        empty_table = LogTailTable(Float32[0.0f0], Float32[0.0f0])
        empty_forward = RaggedArray(Float32[], copy(bundle.forward.offsets))
        empty_reverse = RaggedArray(Float32[], copy(bundle.reverse.offsets))
        return empty_table, StrandPair(empty_forward, empty_reverse)
    end

    flat = Vector{Float32}(undef, n)
    copyto!(flat, 1, bundle.forward.data, 1, n_forward)
    copyto!(flat, n_forward + 1, bundle.reverse.data, 1, n_reverse)
    perm = sortperm(flat; rev=true)

    n_unique = 1
    @inbounds for k in 2:n
        n_unique += flat[perm[k]] != flat[perm[k - 1]]
    end
    unique_scores = Vector{Float32}(undef, n_unique)
    log_tail = Vector{Float32}(undef, n_unique)
    normalized = Vector{Float32}(undef, n)
    group = 1
    k = 1
    @inbounds while k <= n
        score = flat[perm[k]]
        unique_scores[group] = score
        j = k + 1
        while j <= n && flat[perm[j]] == score
            j += 1
        end
        value = Float32(-log10(Float64(j - 1) / Float64(n)))
        log_tail[group] = value
        for p in k:(j - 1)
            normalized[perm[p]] = value
        end
        group += 1
        k = j
    end

    forward_data = @view normalized[1:n_forward]
    reverse_data = @view normalized[(n_forward + 1):n]
    forward = RaggedArray(forward_data, copy(bundle.forward.offsets))
    reverse = RaggedArray(reverse_data, copy(bundle.reverse.offsets))
    return LogTailTable(unique_scores, log_tail), StrandPair(forward, reverse)
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
function transform_scores(table::LogTailTable, scores::RaggedArray{Float32})
    n = length(scores.data)
    n == 0 && return RaggedArray(Float32[], copy(scores.offsets))
    out_data = Vector{Float32}(undef, n)
    @inbounds for i in 1:n
        out_data[i] = lookup_score(table, scores.data[i])
    end
    return RaggedArray(out_data, copy(scores.offsets))
end

"""
    flatten_bundle(bundle::StrandPair{<:RaggedArray{Float32}})

Flatten all valid scores from both strands into a single vector.
Used to fit the normalization table from background scan scores.
"""
function flatten_bundle(bundle::StrandPair{<:RaggedArray{Float32}})
    return vcat(bundle.forward.data, bundle.reverse.data)
end

"""
    normalize_bundle(table::LogTailTable, bundle::StrandPair{<:RaggedArray{Float32}})

Apply the log-tail lookup to both strands of a profile bundle.
"""
function normalize_bundle(table::LogTailTable, bundle::StrandPair{<:RaggedArray{Float32}})
    fwd = transform_scores(table, bundle.forward)
    rev = transform_scores(table, bundle.reverse)
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
