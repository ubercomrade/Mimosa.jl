# Site extraction, PFM reconstruction, and typed selectors.

"""
    SiteSelector

Abstract supertype for site selection strategies.

Concrete selectors:
- [`BestPerSequence`](@ref): one best hit per sequence.
- [`ThresholdHits`](@ref): all hits above a score threshold.
- [`TopFractionHits`](@ref): top fraction of hits by score.
"""
abstract type SiteSelector end

"""
    BestPerSequence

Select the single best-scoring hit per sequence (across the requested strands).
"""
struct BestPerSequence <: SiteSelector end

"""
    ThresholdHits

Select all hits with score ≥ `threshold`.
"""
struct ThresholdHits <: SiteSelector
    threshold::Float32
end

"""
    TopFractionHits

Keep only the top `fraction` of hits by score after collection. At least one
hit is kept. This wraps an underlying selection mode (`best` or `threshold`).
"""
struct TopFractionHits{S<:SiteSelector} <: SiteSelector
    fraction::Float64
    base::S
end

function TopFractionHits(fraction::Real, base::S) where {S<:SiteSelector}
    value = Float64(fraction)
    isfinite(value) && 0.0 < value <= 1.0 ||
        throw(ArgumentError("fraction must be finite and lie in (0, 1]."))
    return TopFractionHits{S}(value, base)
end

function TopFractionHits(fraction::Real)
    return TopFractionHits(fraction, BestPerSequence())
end

"""
    SiteHit

One motif binding site hit.

Fields:
- `seq_index::Int`: 1-based sequence index.
- `start::Int`: 1-based start position in the sequence (inclusive).
- `strand::Int8`: `0` for forward, `1` for reverse.
- `score::Float32`: scan score at this position.
"""
struct SiteHit
    seq_index::Int
    start::Int
    strand::Int8
    score::Float32
end

"""
    SiteCollection

Typed collection of motif hits with parallel arrays.

Fields:
- `seq_indices::Vector{Int}`: 1-based sequence indices.
- `starts::Vector{Int}`: 1-based start positions (inclusive).
- `strands::Vector{Int8}`: `0` for forward, `1` for reverse.
- `scores::Vector{Float32}`: scan scores.
"""
struct SiteCollection
    seq_indices::Vector{Int}
    starts::Vector{Int}
    strands::Vector{Int8}
    scores::Vector{Float32}

    function SiteCollection(
        seq_indices::Vector{Int},
        starts::Vector{Int},
        strands::Vector{Int8},
        scores::Vector{Float32},
    )
        n = length(seq_indices)
        length(starts) == n ||
            throw(ArgumentError("starts length $(length(starts)) != $n."))
        length(strands) == n ||
            throw(ArgumentError("strands length $(length(strands)) != $n."))
        length(scores) == n ||
            throw(ArgumentError("scores length $(length(scores)) != $n."))
        all(>(0), seq_indices) || throw(ArgumentError("sequence indices must be positive."))
        all(>(0), starts) || throw(ArgumentError("site starts must be positive."))
        all(strand -> strand == 0 || strand == 1, strands) ||
            throw(ArgumentError("strands must be 0 (forward) or 1 (reverse)."))
        all(isfinite, scores) || throw(ArgumentError("site scores must be finite."))
        return new(seq_indices, starts, strands, scores)
    end
end

Base.length(coll::SiteCollection) = length(coll.seq_indices)
Base.isempty(coll::SiteCollection) = isempty(coll.seq_indices)

function Base.show(io::IO, coll::SiteCollection)
    return print(io, "SiteCollection($(length(coll)) hits)")
end

function Base.:(==)(a::SiteCollection, b::SiteCollection)
    return a.seq_indices == b.seq_indices &&
           a.starts == b.starts &&
           a.strands == b.strands &&
           a.scores == b.scores
end

"""
    empty_site_collection()

Return an empty [`SiteCollection`](@ref).
"""
empty_site_collection() = SiteCollection(Int[], Int[], Int8[], Float32[])

"""
    SiteHit(coll::SiteCollection, i::Int)

Return the `i`-th hit as a [`SiteHit`](@ref).
"""
function SiteHit(coll::SiteCollection, i::Int)
    return SiteHit(coll.seq_indices[i], coll.starts[i], coll.strands[i], coll.scores[i])
end

# ── Best-per-sequence selection ──────────────────────────────────────────

"""
    _collect_best_hits(bundle::StrandPair{<:RaggedArray{Float32}})

Collect the single best hit per sequence across both strands, matching
Python's `_collect_best_hits`. Returns a [`SiteCollection`](@ref).

Tie-breaking: strictly greater-than comparison, so on equal scores the
forward strand (strand 0) wins, and within a strand the first position wins
(argmax returns the first maximum).
"""
function _collect_best_hits(bundle::StrandPair{<:RaggedArray{Float32}})
    n = nrows(bundle.forward)
    seq_indices = Int[]
    starts = Int[]
    strands = Int8[]
    scores = Float32[]

    for seq_idx in 1:n
        best_score = Float32(-Inf)
        best_start = -1
        best_strand = Int8(0)

        # Forward strand (strand 0)
        fwd_scores = row(bundle.forward, seq_idx)
        if length(fwd_scores) > 0
            fwd_start, fwd_score = _argmax_first(fwd_scores)
            if fwd_score > best_score
                best_score = fwd_score
                best_start = fwd_start
                best_strand = Int8(0)
            end
        end

        # Reverse strand (strand 1)
        rev_scores = row(bundle.reverse, seq_idx)
        if length(rev_scores) > 0
            rev_start, rev_score = _argmax_first(rev_scores)
            if rev_score > best_score
                best_score = rev_score
                best_start = rev_start
                best_strand = Int8(1)
            end
        end

        if best_start < 0 || !isfinite(best_score)
            continue
        end

        push!(seq_indices, seq_idx)
        push!(starts, best_start)
        push!(strands, best_strand)
        push!(scores, best_score)
    end

    return SiteCollection(seq_indices, starts, strands, scores)
end

"""
    _argmax_first(v)

Return `(index, value)` of the first maximum element in `v` (1-based index).
"""
function _argmax_first(v::AbstractVector{T}) where {T}
    best_idx = 1
    best_val = v[1]
    @inbounds for i in 2:length(v)
        if v[i] > best_val
            best_val = v[i]
            best_idx = i
        end
    end
    return (best_idx, best_val)
end

# ── Threshold selection ─────────────────────────────────────────────────

"""
    _collect_threshold_hits(bundle::StrandPair{<:RaggedArray{Float32}}, threshold::Float32)

Collect all hits with score ≥ `threshold` from both strands.
"""
function _collect_threshold_hits(
    bundle::StrandPair{<:RaggedArray{Float32}}, threshold::Float32
)
    n = nrows(bundle.forward)
    seq_indices = Int[]
    starts = Int[]
    strands = Int8[]
    scores = Float32[]

    for seq_idx in 1:n
        # Forward strand
        fwd_scores = row(bundle.forward, seq_idx)
        for pos in 1:length(fwd_scores)
            if fwd_scores[pos] >= threshold
                push!(seq_indices, seq_idx)
                push!(starts, pos)
                push!(strands, Int8(0))
                push!(scores, fwd_scores[pos])
            end
        end

        # Reverse strand
        rev_scores = row(bundle.reverse, seq_idx)
        for pos in 1:length(rev_scores)
            if rev_scores[pos] >= threshold
                push!(seq_indices, seq_idx)
                push!(starts, pos)
                push!(strands, Int8(1))
                push!(scores, rev_scores[pos])
            end
        end
    end

    return SiteCollection(seq_indices, starts, strands, scores)
end

# ── Best-strand threshold selection ──────────────────────────────────────

"""
    _collect_best_strand_threshold_hits(bundle::StrandPair{<:RaggedArray{Float32}}, threshold::Float32)

Collect above-threshold hits after collapsing both strands by per-position maximum.
The strand is chosen as forward if `fwd >= rev` at that position, matching Python's
`_collect_best_strand_threshold_hits`.
"""
function _collect_best_strand_threshold_hits(
    bundle::StrandPair{<:RaggedArray{Float32}}, threshold::Float32
)
    n = nrows(bundle.forward)
    seq_indices = Int[]
    starts = Int[]
    strands = Int8[]
    scores = Float32[]

    for seq_idx in 1:n
        fwd_scores = row(bundle.forward, seq_idx)
        rev_scores = row(bundle.reverse, seq_idx)
        if isempty(fwd_scores)
            continue
        end

        for pos in 1:length(fwd_scores)
            best = max(fwd_scores[pos], rev_scores[pos])
            if best >= threshold
                push!(seq_indices, seq_idx)
                push!(starts, pos)
                push!(strands, Int8(fwd_scores[pos] >= rev_scores[pos] ? 0 : 1))
                push!(scores, best)
            end
        end
    end

    return SiteCollection(seq_indices, starts, strands, scores)
end

# ── Sort hits ────────────────────────────────────────────────────────────

"""
    sort_hits!(coll::SiteCollection)

Sort hits in-place by (seq_index ascending, score descending, start ascending,
strand ascending), matching Python's `_sort_hit_arrays` via `np.lexsort`.
"""
function sort_hits!(coll::SiteCollection)
    n = length(coll)
    n <= 1 && return coll

    # Build sort keys: primary = seq_index asc, secondary = -score asc (score desc),
    # tertiary = start asc, quaternary = strand asc
    # Julia's sortperm with by/lt or sortslices is tricky for multi-key.
    # Use a simple approach: sort by tuple.
    perm = sortperm(
        1:n; by=i -> (coll.seq_indices[i], -coll.scores[i], coll.starts[i], coll.strands[i])
    )

    coll.seq_indices[:] = coll.seq_indices[perm]
    coll.starts[:] = coll.starts[perm]
    coll.strands[:] = coll.strands[perm]
    coll.scores[:] = coll.scores[perm]
    return coll
end

# ── Top-fraction selection ───────────────────────────────────────────────

"""
    select_top_fraction(coll::SiteCollection, fraction::Float64)

Keep only the top `fraction` of hits by score. At least one hit is kept.
Matches Python's `_select_top_hit_arrays`.
"""
function select_top_fraction(coll::SiteCollection, fraction::Float64)
    n = length(coll)
    n == 0 && return coll
    isfinite(fraction) && 0.0 < fraction <= 1.0 ||
        throw(ArgumentError("fraction must be finite and lie in (0, 1]."))

    n_keep = max(1, floor(Int, n * fraction))
    n_keep >= n && return coll

    # Find the n_keep highest-scoring indices
    perm = sortperm(coll.scores; rev=true)
    keep = perm[1:n_keep]
    # Re-sort kept hits by score descending
    keep = keep[sortperm(coll.scores[keep]; rev=true)]

    return SiteCollection(
        coll.seq_indices[keep], coll.starts[keep], coll.strands[keep], coll.scores[keep]
    )
end

# ── Site extraction ──────────────────────────────────────────────────────

"""
    extract_site_matrix(batch::EncodedSequenceBatch, coll::SiteCollection, motif_width::Int;
                         site_offset::Int=0)

Extract numeric site windows for all hits. For reverse-strand hits, the site
is reverse-complemented so that it is in canonical forward motif orientation.

`site_offset` is the number of bases between the scan position and the actual
motif start (e.g. `context_length` for BaMM/Dimont/Slim). For PWM/PFM/SiteGA
this is 0.

Returns a `Matrix{UInt8}` of shape `(motif_width, n_hits)` where each column
is one site window.
"""
function extract_site_matrix(
    batch::EncodedSequenceBatch, coll::SiteCollection, motif_width::Int; site_offset::Int=0
)
    n_hits = length(coll)
    motif_width > 0 || throw(ArgumentError("motif_width must be positive."))
    site_offset >= 0 || throw(ArgumentError("site_offset must be non-negative."))
    nsequences(batch) >= 0 || throw(InvariantError("invalid sequence batch."))
    all(1 <= i <= nsequences(batch) for i in coll.seq_indices) ||
        throw(ArgumentError("site sequence indices are outside the batch."))
    all(strand -> strand == 0 || strand == 1, coll.strands) ||
        throw(ArgumentError("site strands must be 0 or 1."))
    sites = Matrix{UInt8}(undef, motif_width, n_hits)

    for h in 1:n_hits
        seq = sequence(batch, coll.seq_indices[h])
        # The motif starts at scan_position + site_offset (1-based)
        coll.starts[h] <= typemax(Int) - site_offset ||
            throw(ArgumentError("site start plus offset overflows Int."))
        start = coll.starts[h] + site_offset
        # Validate that the site window fits within the sequence.
        if start < 1 || start > length(seq) || motif_width > length(seq) - start + 1
            throw(
                InvariantError(
                    "site $h: window [start=$start, width=$motif_width] " *
                    "exceeds sequence length $(length(seq)).",
                ),
            )
        end
        # Extract window: seq[start : start + motif_width - 1] (1-based)
        # Invariant: seq codes in 0..N_CODE (guaranteed by EncodedSequenceBatch).
        # @inbounds is safe: start >= 1 and start+motif_width-1 <= length(seq).
        @inbounds for p in 1:motif_width
            sites[p, h] = seq[start + p - 1]
        end

        # Reverse complement for minus strand
        if coll.strands[h] == 1
            _reverse_complement_site!(view(sites, :, h), motif_width)
        end
    end

    return sites
end

"""
    _reverse_complement_site!(site::AbstractVector{UInt8}, len::Int)

Reverse-complement a site window in-place. N (0x04) stays as N.
"""
function _reverse_complement_site!(site::AbstractVector{UInt8}, len::Int)
    @inbounds for i in 1:(len ÷ 2)
        j = len - i + 1
        a = site[i]
        b = site[j]
        site[i] = b == N_CODE ? N_CODE : 0x03 - b
        site[j] = a == N_CODE ? N_CODE : 0x03 - a
    end
    if isodd(len)
        mid = (len + 1) ÷ 2
        @inbounds site[mid] = site[mid] == N_CODE ? N_CODE : 0x03 - site[mid]
    end
    return site
end

# ── Site strings ─────────────────────────────────────────────────────────

const _SEQ_DECODER = UInt8['A', 'C', 'G', 'T', 'N']

"""
    site_strings(sites::Matrix{UInt8})

Convert numeric site windows to DNA strings. Each column of `sites` is one
site window. Returns a `Vector{String}`.
"""
function site_strings(sites::Matrix{UInt8})
    n_hits = size(sites, 2)
    motif_width = size(sites, 1)
    result = Vector{String}(undef, n_hits)
    for h in 1:n_hits
        buf = Vector{UInt8}(undef, motif_width)
        @inbounds for p in 1:motif_width
            code = sites[p, h]
            idx = code > N_CODE ? N_CODE + 1 : code + 1
            buf[p] = _SEQ_DECODER[idx]
        end
        result[h] = String(buf)
    end
    return result
end

# ── PCM and PFM reconstruction ───────────────────────────────────────────

"""
    build_pcm(sites::Matrix{UInt8}, motif_width::Int)

Build a Position Count Matrix from numeric site windows.

Only valid bases (A=0, C=1, G=2, T=3) are counted; N (4) is skipped.
Returns a `Matrix{Float32}` of shape `(4, motif_width)`.
"""
function build_pcm(sites::Matrix{UInt8}, motif_width::Int)
    motif_width > 0 || throw(ArgumentError("motif_width must be positive."))
    size(sites, 1) == motif_width ||
        throw(ArgumentError("sites row count must equal motif_width."))
    all(code -> code <= N_CODE, sites) ||
        throw(ArgumentError("sites contain invalid DNA codes."))
    pcm = zeros(Float32, 4, motif_width)
    n_hits = size(sites, 2)
    @inbounds for h in 1:n_hits
        for p in 1:motif_width
            code = sites[p, h]
            if code < N_CODE
                pcm[code + 1, p] += 1.0f0
            end
        end
    end
    return pcm
end

"""
    pcm_to_pfm(pcm; pseudocount=0.25)

Convert a Position Count Matrix to a Position Frequency Matrix.

`pcm` axes: `(base, position)` with `base ∈ 1:4`.
"""
function pcm_to_pfm(
    pcm::AbstractMatrix{T}; pseudocount::AbstractFloat=0.25f0
) where {T<:AbstractFloat}
    if size(pcm, 1) != NUCLEOTIDE_CARDINALITY
        throw(ModelDimensionError("PCM must have 4 rows, got $(size(pcm, 1))."))
    end
    n_sites = sum(pcm; dims=1)
    pc = T(pseudocount)
    denom = n_sites .+ T(NUCLEOTIDE_CARDINALITY) * pc
    return (pcm .+ pc) ./ denom
end

"""
    _empty_ragged_like(rag::RaggedArray{T}) where T

Return a `RaggedArray{T}` with the same number of rows but all empty (zero-length).
"""
function _empty_ragged_like(rag::RaggedArray{T}) where {T}
    n = nrows(rag)
    return RaggedArray(Vector{T}(), fill(1, n + 1))
end

# ── Generic AbstractMotifModel support (Extensibility API Plan §7.2) ─────────
#
# The hit collection logic (best, threshold, top-fraction) is identical across
# all model families — it operates on StrandPair score bundles. The only
# differences are:
#   1. How the bundle is produced (dispatch via scan(model, batch; ...)).
#   2. The site start offset (left_context for any model with context).
# Both are resolved via the public geometry accessors.

"""
    _scan_bundle_for_sites(model::AbstractMotifModel, batch, strands)

Scan any `AbstractMotifModel` (custom or higher-order built-in) to produce the
required strand bundle.
"""
function _scan_bundle_for_sites(
    model::AbstractMotifModel,
    batch::EncodedSequenceBatch,
    strands::StrandPolicy,
    execution::ExecutionPolicy=SerialExecution(),
)
    if strands isa ForwardOnly
        fwd = _scan_model_batch(model, batch; strands=ForwardOnly(), execution=execution)
        rev = _empty_ragged_like(fwd)
        return StrandPair(fwd, rev)
    elseif strands isa ReverseOnly
        rev = _scan_model_batch(model, batch; strands=ReverseOnly(), execution=execution)
        fwd = _empty_ragged_like(rev)
        return StrandPair(fwd, rev)
    else
        result = _scan_model_batch(model, batch; strands=BothStrands(), execution=execution)
        return result
    end
end

"""
    _collect_hits(model::AbstractMotifModel, batch, selector; strands)

Collect hits from any model scan. The hit collection logic is shared with PWM
since it operates on StrandPair bundles.
"""
function _collect_hits(
    model::AbstractMotifModel,
    batch::EncodedSequenceBatch,
    selector::SiteSelector;
    strands::StrandPolicy=BothStrands(),
    execution::ExecutionPolicy=SerialExecution(),
)
    bundle = _scan_bundle_for_sites(model, batch, strands, execution)
    return _collect_hits_from_bundle(selector, bundle, strands)
end

# Shared dispatch from bundle to hits (no model dependency).
function _collect_hits_from_bundle(
    selector::SiteSelector,
    bundle::StrandPair{<:RaggedArray{Float32}},
    strands::StrandPolicy,
)
    if selector isa BestPerSequence
        return _collect_best_hits(bundle)
    elseif selector isa ThresholdHits
        if strands isa BestStrand
            return _collect_best_strand_threshold_hits(bundle, selector.threshold)
        else
            return _collect_threshold_hits(bundle, selector.threshold)
        end
    elseif selector isa TopFractionHits
        base_coll = _collect_hits_from_bundle(selector.base, bundle, strands)
        return select_top_fraction(base_coll, selector.fraction)
    else
        throw(ArgumentError("unsupported selector: $selector."))
    end
end

function _resolve_selector(mode::Symbol; score_threshold::Union{Nothing,Float32}=nothing)
    if mode == :best
        return BestPerSequence()
    elseif mode == :threshold
        score_threshold === nothing &&
            throw(ArgumentError("score_threshold is required for mode=:threshold."))
        return ThresholdHits(score_threshold)
    else
        throw(ArgumentError("mode must be :best or :threshold, got :$mode."))
    end
end

"""
    selectsites(model::AbstractMotifModel, batch, selector;
                strands=BothStrands(), execution=SerialExecution())

Extract motif binding sites from any motif model. The `start` field in the
returned [`SiteCollection`](@ref) is the scan position; the actual motif start
is `start + site_start_offset(model)` (= `left_context(model)`).

Returns a sorted [`SiteCollection`](@ref).
"""
function selectsites(
    model::AbstractMotifModel,
    batch::EncodedSequenceBatch,
    selector::SiteSelector;
    strands::StrandPolicy=BothStrands(),
    execution::ExecutionPolicy=SerialExecution(),
)
    validate_model(model; capability=:sites)
    coll = _collect_hits(model, batch, selector; strands=strands, execution=execution)
    sort_hits!(coll)
    return coll
end

"""
    reconstruct_pfm(model::AbstractMotifModel, batch, selector;
                    pseudocount=0.25f0, strands=BothStrands(),
                    execution=SerialExecution())

Reconstruct a PFM from binding sites extracted from any motif model. The site
window accounts for `site_start_offset(model)` (= `left_context(model)`) to
extract only the motif-length window, excluding context bases.

Returns a `Matrix{Float32}` of shape `(4, motif_length)`.
"""
function reconstruct_pfm(
    model::AbstractMotifModel,
    batch::EncodedSequenceBatch,
    selector::SiteSelector;
    pseudocount::Float32=0.25f0,
    strands::StrandPolicy=BothStrands(),
    execution::ExecutionPolicy=SerialExecution(),
)
    validate_model(model; capability=:sites)
    coll = _collect_hits(model, batch, selector; strands=strands, execution=execution)

    if isempty(coll)
        throw(ArgumentError("No sites found for PFM reconstruction."))
    end

    motif_width = motif_length(model)
    offset = site_start_offset(model)
    sites = extract_site_matrix(batch, coll, motif_width; site_offset=offset)
    pcm = build_pcm(sites, motif_width)
    return pcm_to_pfm(pcm; pseudocount=pseudocount)
end

"""
    reconstruct_pfm(model::AbstractMotifModel, batch;
                    mode=:best, pseudocount=0.25f0, strands=BothStrands(),
                    score_threshold=nothing, top_fraction=nothing)

Convenience method for PFM reconstruction from any motif model with
keyword-based selection mode.
"""
function reconstruct_pfm(
    model::AbstractMotifModel,
    batch::EncodedSequenceBatch;
    mode::Symbol=:best,
    pseudocount::Float32=0.25f0,
    strands::StrandPolicy=BothStrands(),
    score_threshold::Union{Nothing,Float32}=nothing,
    top_fraction::Union{Nothing,Float64}=nothing,
    execution::ExecutionPolicy=SerialExecution(),
)
    selector = _resolve_selector(mode; score_threshold=score_threshold)
    if top_fraction !== nothing
        selector = TopFractionHits(top_fraction, selector)
    end

    return reconstruct_pfm(
        model,
        batch,
        selector;
        pseudocount=pseudocount,
        strands=strands,
        execution=execution,
    )
end
