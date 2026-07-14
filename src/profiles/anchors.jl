# Anchor collection for profile comparison: best-per-row or threshold-based.

"""
    AnchorCSR

Compressed sparse row representation of anchor positions.

Fields:
- `positions::Vector{Int}`: anchor positions sorted by row (stable order).
- `offsets::Vector{Int}`: row offsets (length `n_rows + 1`), so
  `positions[offsets[i]:offsets[i+1]-1]` are the anchors for row `i`.
"""
struct AnchorCSR
    positions::Vector{Int}
    offsets::Vector{Int}

    function AnchorCSR(positions::Vector{Int}, offsets::Vector{Int})
        isempty(offsets) && throw(ArgumentError("anchor offsets must not be empty."))
        offsets[1] == 1 || throw(ArgumentError("anchor offsets must start at 1."))
        offsets[end] == length(positions) + 1 ||
            throw(ArgumentError("anchor offsets must end at positions length + 1."))
        all(diff(offsets) .>= 0) ||
            throw(ArgumentError("anchor offsets must be nondecreasing."))
        return new(positions, offsets)
    end
end

Base.isempty(csr::AnchorCSR) = isempty(csr.positions)

"""
    build_anchor_csr(rows::Vector{Int}, positions::Vector{Int}, n_rows::Int)

Build an [`AnchorCSR`](@ref) from flat row/position arrays using a stable
counting sort. Matches Python's `build_anchor_csr`.
"""
function build_anchor_csr(rows::Vector{Int}, positions::Vector{Int}, n_rows::Int)
    n_rows >= 0 || throw(ArgumentError("n_rows must be non-negative."))
    length(rows) == length(positions) ||
        throw(ArgumentError("rows and positions must have equal lengths."))
    all(1 <= row <= n_rows for row in rows) ||
        throw(ArgumentError("anchor rows must be within 1:n_rows."))
    all(position > 0 for position in positions) ||
        throw(ArgumentError("anchor positions must be positive."))

    n = length(rows)
    n == 0 && return AnchorCSR(Int[], ones(Int, n_rows + 1))

    # Count per row
    counts = zeros(Int, n_rows)
    @inbounds for r in rows
        counts[r] += 1
    end

    # Build offsets (1-based)
    offsets = Vector{Int}(undef, n_rows + 1)
    offsets[1] = 1
    @inbounds for i in 1:n_rows
        offsets[i + 1] = offsets[i] + counts[i]
    end

    # Stable counting sort: place positions preserving input order within each row
    sorted_positions = Vector{Int}(undef, n)
    next_slot = copy(offsets)
    @inbounds for i in 1:n
        sorted_positions[next_slot[rows[i]]] = positions[i]
        next_slot[rows[i]] += 1
    end

    return AnchorCSR(sorted_positions, offsets)
end

"""
    collect_best_anchors(scores::RaggedArray{Float32})

Collect the best (highest-scoring) position per non-empty row.
Returns `(rows, positions)` as flat vectors. First position wins on ties
(matching Python's `_collect_best_anchor_positions_numba`).
"""
function collect_best_anchors(scores::RaggedArray{Float32})
    n = nrows(scores)
    rows = Int[]
    positions = Int[]
    for i in 1:n
        len = rowlength(scores, i)
        len == 0 && continue
        r = row(scores, i)
        best_pos = 1
        best_score = r[1]
        @inbounds for j in 2:len
            if r[j] > best_score
                best_score = r[j]
                best_pos = j
            end
        end
        push!(rows, i)
        push!(positions, best_pos)
    end
    return (rows, positions)
end

"""
    collect_threshold_anchors(scores::RaggedArray{Float32}, threshold::Float32)

Collect all positions where `score >= threshold`. Returns `(rows, positions)`.
"""
function collect_threshold_anchors(scores::RaggedArray{Float32}, threshold::Float32)
    n = nrows(scores)
    rows = Int[]
    positions = Int[]
    for i in 1:n
        r = row(scores, i)
        for j in eachindex(r)
            if r[j] >= threshold
                push!(rows, i)
                push!(positions, j)
            end
        end
    end
    return (rows, positions)
end

"""
    collect_anchors(scores::RaggedArray{Float32}, threshold::Float32)

Dispatch to best or threshold anchor collection based on the threshold value.
A threshold `<= 0` means "use best anchors" (matching Python's
`score_threshold = None if min_logfpr is None or min_logfpr <= 0`).
"""
function collect_anchors(scores::RaggedArray{Float32}, threshold::Float32)
    if threshold > 0
        return collect_threshold_anchors(scores, threshold)
    else
        return collect_best_anchors(scores)
    end
end
