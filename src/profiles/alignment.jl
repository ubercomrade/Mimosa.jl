# Profile window alignment: shift-based comparison with anchor realignment.

"""
    _window_fits(pos::Int, len::Int, radius::Int)

Check if a window of radius `radius` centered at 1-based position `pos` fits
within `[1, len]`. Equivalent to Python's `_window_fits` (0-based: `pos - radius >= 0 && pos + radius < length`).
"""
function _window_fits(pos::Int, len::Int, radius::Int)
    return pos - radius >= 1 && pos + radius <= len
end

"""
    _realign_query_position(r::AbstractVector{Float32}, expected::Int, radius::Int)

Find the best-scoring position in `[max(1, expected-radius), min(len, expected+radius)]`.
Returns 0 if the range is empty. Matches Python's `_realign_query_position`.
"""
function _realign_query_position(r::AbstractVector{Float32}, expected::Int, radius::Int)
    len = length(r)
    left = max(1, expected - radius)
    right = min(len, expected + radius)
    left > right && return 0
    best_pos = left
    best_score = r[left]
    @inbounds for pos in (left + 1):right
        if r[pos] > best_score
            best_score = r[pos]
            best_pos = pos
        end
    end
    return best_pos
end

mutable struct CandidateScratch
    candidates::Vector{Int}
    seen_epoch::Vector{UInt32}
    epoch::UInt32
end

CandidateScratch(max_len::Int) = CandidateScratch(Int[], zeros(UInt32, max_len), 0)

function _collect_row_candidates!(
    scratch::CandidateScratch,
    r1::AbstractVector{Float32},
    len1::Int,
    len2::Int,
    query_csr::AnchorCSR,
    target_csr::AnchorCSR,
    row::Int,
    shift::Int,
    window_radius::Int,
    realign_window::Int,
)
    empty!(scratch.candidates)
    if length(scratch.seen_epoch) < len1
        resize!(scratch.seen_epoch, len1)
        fill!(scratch.seen_epoch, 0)
        scratch.epoch = 0
    end
    if scratch.epoch == typemax(UInt32)
        fill!(scratch.seen_epoch, 0)
        scratch.epoch = UInt32(1)
    else
        scratch.epoch += UInt32(1)
    end
    epoch = scratch.epoch

    # Query anchors
    @inbounds for idx in query_csr.offsets[row]:(query_csr.offsets[row + 1] - 1)
        pos1 = query_csr.positions[idx]
        pos2 = pos1 + shift
        if _window_fits(pos1, len1, window_radius) &&
            _window_fits(pos2, len2, window_radius)
            if scratch.seen_epoch[pos1] != epoch
                scratch.seen_epoch[pos1] = epoch
                push!(scratch.candidates, pos1)
            end
        end
    end

    # Target anchors (realigned to query coordinates)
    @inbounds for idx in target_csr.offsets[row]:(target_csr.offsets[row + 1] - 1)
        expected_pos1 = target_csr.positions[idx] - shift
        pos1 = _realign_query_position(r1, expected_pos1, realign_window)
        pos1 == 0 && continue
        pos2 = pos1 + shift
        if _window_fits(pos1, len1, window_radius) &&
            _window_fits(pos2, len2, window_radius)
            if scratch.seen_epoch[pos1] != epoch
                scratch.seen_epoch[pos1] = epoch
                push!(scratch.candidates, pos1)
            end
        end
    end

    return scratch.candidates
end

# Accumulate pooled overlap sums (for CO and Dice metrics).
function _accumulate_pooled(
    r1::AbstractVector{Float32},
    r2::AbstractVector{Float32},
    candidates::Vector{Int},
    shift::Int,
    window_radius::Int,
)
    sum1 = 0.0
    sum2 = 0.0
    intersection = 0.0
    @inbounds for pos1 in candidates
        pos2 = pos1 + shift
        for offset in (-window_radius):window_radius
            v1 = Float64(r1[pos1 + offset])
            v2 = Float64(r2[pos2 + offset])
            sum1 += v1
            sum2 += v2
            intersection += min(v1, v2)
        end
    end
    return (sum1, sum2, intersection)
end

function _accumulate_pooled(
    r1::AbstractVector{Float32},
    r2::AbstractVector{Float32},
    pos1::Int,
    shift::Int,
    window_radius::Int,
)
    pos2 = pos1 + shift
    sum1 = 0.0
    sum2 = 0.0
    intersection = 0.0
    @inbounds for offset in (-window_radius):window_radius
        v1 = Float64(r1[pos1 + offset])
        v2 = Float64(r2[pos2 + offset])
        sum1 += v1
        sum2 += v2
        intersection += min(v1, v2)
    end
    return (sum1, sum2, intersection)
end

# Accumulate rowwise overlap-based values (for CO_rowwise and Dice_rowwise).
function _accumulate_rowwise_overlap(
    r1::AbstractVector{Float32},
    r2::AbstractVector{Float32},
    candidates::Vector{Int},
    shift::Int,
    window_radius::Int,
    use_dice::Bool,
)
    score_sum = 0.0
    finite_count = 0
    @inbounds for pos1 in candidates
        pos2 = pos1 + shift
        sum1 = 0.0
        sum2 = 0.0
        intersection = 0.0
        for offset in (-window_radius):window_radius
            v1 = Float64(r1[pos1 + offset])
            v2 = Float64(r2[pos2 + offset])
            sum1 += v1
            sum2 += v2
            intersection += min(v1, v2)
        end
        denom = use_dice ? sum1 + sum2 : min(sum1, sum2)
        if denom > PROFILE_EPS
            val = use_dice ? 2.0 * intersection / denom : intersection / denom
            score_sum += val
            finite_count += 1
        end
    end
    return (score_sum, finite_count)
end

function _accumulate_rowwise_overlap(
    r1::AbstractVector{Float32},
    r2::AbstractVector{Float32},
    pos1::Int,
    shift::Int,
    window_radius::Int,
    use_dice::Bool,
)
    pos2 = pos1 + shift
    sum1 = 0.0
    sum2 = 0.0
    intersection = 0.0
    @inbounds for offset in (-window_radius):window_radius
        v1 = Float64(r1[pos1 + offset])
        v2 = Float64(r2[pos2 + offset])
        sum1 += v1
        sum2 += v2
        intersection += min(v1, v2)
    end
    denom = use_dice ? sum1 + sum2 : min(sum1, sum2)
    if denom > PROFILE_EPS
        return (use_dice ? 2.0 * intersection / denom : intersection / denom, 1)
    end
    return (0.0, 0)
end

# Accumulate rowwise cosine values.
function _accumulate_cosine(
    r1::AbstractVector{Float32},
    r2::AbstractVector{Float32},
    candidates::Vector{Int},
    shift::Int,
    window_radius::Int,
)
    score_sum = 0.0
    finite_count = 0
    @inbounds for pos1 in candidates
        pos2 = pos1 + shift
        dot = 0.0
        norm1 = 0.0
        norm2 = 0.0
        for offset in (-window_radius):window_radius
            v1 = Float64(r1[pos1 + offset])
            v2 = Float64(r2[pos2 + offset])
            dot += v1 * v2
            norm1 += v1 * v1
            norm2 += v2 * v2
        end
        denom = sqrt(norm1) * sqrt(norm2)
        if denom > PROFILE_EPS
            score_sum += dot / denom
            finite_count += 1
        end
    end
    return (score_sum, finite_count)
end

function _accumulate_cosine(
    r1::AbstractVector{Float32},
    r2::AbstractVector{Float32},
    pos1::Int,
    shift::Int,
    window_radius::Int,
)
    pos2 = pos1 + shift
    dot = 0.0
    norm1 = 0.0
    norm2 = 0.0
    @inbounds for offset in (-window_radius):window_radius
        v1 = Float64(r1[pos1 + offset])
        v2 = Float64(r2[pos2 + offset])
        dot += v1 * v2
        norm1 += v1 * v1
        norm2 += v2 * v2
    end
    denom = sqrt(norm1) * sqrt(norm2)
    return denom > PROFILE_EPS ? (dot / denom, 1) : (0.0, 0)
end

"""
    score_shift(scores1, scores2, query_csr, target_csr, shift, window_radius, realign_window, metric)

Evaluate one shift across all rows and reduce row-level partials into a
`(score::Float32, n_sites::Int)` tuple. Matches Python's `score_shift` serial path.

- `scores1`: query strand scores (RaggedArray, 1-based positions per row).
- `scores2`: target strand scores.
- `query_csr`/`target_csr`: anchor CSR structures for this orientation pair.
- `shift`: target position = query position + shift.
- `window_radius`: half-window size (0 = single position).
- `realign_window`: realignment search radius for target anchors.
- `metric`: typed profile metric.
"""
function score_shift(
    scores1::RaggedArray{Float32},
    scores2::RaggedArray{Float32},
    query_csr::AnchorCSR,
    target_csr::AnchorCSR,
    shift::Int,
    window_radius::Int,
    realign_window::Int,
    metric::AbstractProfileMetric,
)
    max_len = maximum((rowlength(scores1, r) for r in 1:nrows(scores1)); init=0)
    scratch = CandidateScratch(max_len)
    return _score_shift!(
        scratch,
        scores1,
        scores2,
        query_csr,
        target_csr,
        shift,
        window_radius,
        realign_window,
        metric,
    )
end

function _score_shift!(
    scratch::CandidateScratch,
    scores1::RaggedArray{Float32},
    scores2::RaggedArray{Float32},
    query_csr::AnchorCSR,
    target_csr::AnchorCSR,
    shift::Int,
    window_radius::Int,
    realign_window::Int,
    metric::AbstractProfileMetric,
)
    n = nrows(scores1)
    total_sum1 = 0.0
    total_sum2 = 0.0
    total_intersection = 0.0
    total_row_score = 0.0
    total_finite = 0
    total_sites = 0

    for r in 1:n
        len1 = rowlength(scores1, r)
        len2 = rowlength(scores2, r)

        r1 = row(scores1, r)
        r2 = row(scores2, r)

        candidates = _collect_row_candidates!(
            scratch,
            r1,
            len1,
            len2,
            query_csr,
            target_csr,
            r,
            shift,
            window_radius,
            realign_window,
        )
        count = length(candidates)
        total_sites += count

        if count == 0
            continue
        end

        if is_pooled(metric)
            s1, s2, inter = _accumulate_pooled(r1, r2, candidates, shift, window_radius)
            total_sum1 += s1
            total_sum2 += s2
            total_intersection += inter
        elseif metric isa CosineSimilarityProfile
            s_sum, f_count = _accumulate_cosine(r1, r2, candidates, shift, window_radius)
            total_row_score += s_sum
            total_finite += f_count
        else
            use_dice = is_dice_metric(metric)
            s_sum, f_count = _accumulate_rowwise_overlap(
                r1, r2, candidates, shift, window_radius, use_dice
            )
            total_row_score += s_sum
            total_finite += f_count
        end
    end

    if total_sites == 0
        return (0.0f0, 0)
    end

    if is_pooled(metric)
        if is_dice_metric(metric)
            denom = total_sum1 + total_sum2
            score = denom > PROFILE_EPS ? 2.0 * total_intersection / denom : 0.0
        else
            denom = min(total_sum1, total_sum2)
            score = denom > PROFILE_EPS ? total_intersection / denom : 0.0
        end
    else
        score = total_finite == 0 ? 0.0 : total_row_score / total_finite
    end

    return (Float32(score), total_sites)
end

function _score_shift_best!(
    scores1::RaggedArray{Float32},
    scores2::RaggedArray{Float32},
    query_csr::AnchorCSR,
    target_csr::AnchorCSR,
    shift::Int,
    window_radius::Int,
    realign_window::Int,
    metric::AbstractProfileMetric,
)
    total_sum1 = 0.0
    total_sum2 = 0.0
    total_intersection = 0.0
    total_row_score = 0.0
    total_finite = 0
    total_sites = 0

    for row_index in 1:nrows(scores1)
        len1 = rowlength(scores1, row_index)
        len2 = rowlength(scores2, row_index)
        r1 = row(scores1, row_index)
        r2 = row(scores2, row_index)
        query_pos = 0
        target_pos = 0

        query_start = query_csr.offsets[row_index]
        query_stop = query_csr.offsets[row_index + 1] - 1
        if query_start <= query_stop
            candidate = query_csr.positions[query_start]
            if _window_fits(candidate, len1, window_radius) &&
                _window_fits(candidate + shift, len2, window_radius)
                query_pos = candidate
            end
        end

        target_start = target_csr.offsets[row_index]
        target_stop = target_csr.offsets[row_index + 1] - 1
        if target_start <= target_stop
            expected = target_csr.positions[target_start] - shift
            candidate = _realign_query_position(r1, expected, realign_window)
            if candidate != 0 &&
                _window_fits(candidate, len1, window_radius) &&
                _window_fits(candidate + shift, len2, window_radius)
                target_pos = candidate
            end
        end

        if query_pos != 0
            total_sites += 1
            if is_pooled(metric)
                s1, s2, inter = _accumulate_pooled(r1, r2, query_pos, shift, window_radius)
                total_sum1 += s1
                total_sum2 += s2
                total_intersection += inter
            elseif metric isa CosineSimilarityProfile
                s_sum, f_count = _accumulate_cosine(r1, r2, query_pos, shift, window_radius)
                total_row_score += s_sum
                total_finite += f_count
            else
                s_sum, f_count = _accumulate_rowwise_overlap(
                    r1, r2, query_pos, shift, window_radius, is_dice_metric(metric)
                )
                total_row_score += s_sum
                total_finite += f_count
            end
        end

        if target_pos != 0 && target_pos != query_pos
            total_sites += 1
            if is_pooled(metric)
                s1, s2, inter = _accumulate_pooled(r1, r2, target_pos, shift, window_radius)
                total_sum1 += s1
                total_sum2 += s2
                total_intersection += inter
            elseif metric isa CosineSimilarityProfile
                s_sum, f_count = _accumulate_cosine(
                    r1, r2, target_pos, shift, window_radius
                )
                total_row_score += s_sum
                total_finite += f_count
            else
                s_sum, f_count = _accumulate_rowwise_overlap(
                    r1, r2, target_pos, shift, window_radius, is_dice_metric(metric)
                )
                total_row_score += s_sum
                total_finite += f_count
            end
        end
    end

    total_sites == 0 && return (0.0f0, 0)
    if is_pooled(metric)
        if is_dice_metric(metric)
            denom = total_sum1 + total_sum2
            score = denom > PROFILE_EPS ? 2.0 * total_intersection / denom : 0.0
        else
            denom = min(total_sum1, total_sum2)
            score = denom > PROFILE_EPS ? total_intersection / denom : 0.0
        end
    else
        score = total_finite == 0 ? 0.0 : total_row_score / total_finite
    end
    return (Float32(score), total_sites)
end

# ── Orientation scoring ──────────────────────────────────────────────────

# Profile orientation pairs: (label, query_strand, target_strand).
# Strand indices: 1=forward, 2=reverse (1-based Julia convention).
const PROFILE_ORIENTATION_PAIRS = (("++", 1, 1), ("+-", 1, 2), ("-+", 2, 1), ("--", 2, 2))
const PROFILE_ORIENTATION_RANK = Dict("++" => 0, "+-" => 1, "-+" => 2, "--" => 3)

function _orientation_pairs(query_bundle, target_bundle)
    qs = query_bundle.forward === query_bundle.reverse
    ts = target_bundle.forward === target_bundle.reverse
    qs && ts && return (("++", 1, 1),)
    qs && return (("++", 1, 1), ("+-", 1, 2))
    ts && return (("++", 1, 1), ("-+", 2, 1))
    return PROFILE_ORIENTATION_PAIRS
end

"""
    _score_orientation_pair(query_bundle, target_bundle, query_anchors, target_anchors, search_range, window_radius, realign_window, metric)

Score one orientation pair across all shifts in `[-search_range, search_range]`.
Returns `(score, shift, n_sites, orientation_label)`.

Tie-breaking (matching Python):
1. Higher score wins.
2. On equal score: more n_sites wins.
3. On equal n_sites: smaller |shift| wins.
4. First in iteration order wins on complete tie.
"""
function _score_orientation_pair_csr(
    query_bundle::StrandPair{<:RaggedArray{Float32}},
    target_bundle::StrandPair{<:RaggedArray{Float32}},
    query_anchors::AnchorCSR,
    target_anchors::AnchorCSR,
    query_strand::Int,
    target_strand::Int,
    orientation_label::String,
    search_range::Int,
    window_radius::Int,
    realign_window::Int,
    metric::AbstractProfileMetric,
)
    query_scores = query_strand == 1 ? query_bundle.forward : query_bundle.reverse
    target_scores = target_strand == 1 ? target_bundle.forward : target_bundle.reverse

    best_score = 0.0f0
    best_shift = 0
    best_n_sites = 0
    max_len = maximum((rowlength(query_scores, r) for r in 1:nrows(query_scores)); init=0)
    scratch = CandidateScratch(max_len)

    for shift in (-search_range):search_range
        score, n_sites = _score_shift!(
            scratch,
            query_scores,
            target_scores,
            query_anchors,
            target_anchors,
            shift,
            window_radius,
            realign_window,
            metric,
        )
        if Float64(score) > Float64(best_score) || (
            Float64(score) == Float64(best_score) && (
                n_sites > best_n_sites ||
                (n_sites == best_n_sites && abs(shift) < abs(best_shift))
            )
        )
            best_score = score
            best_shift = shift
            best_n_sites = n_sites
        end
    end

    return (best_score, best_shift, best_n_sites, orientation_label)
end

function _score_orientation_pair(
    query_bundle::StrandPair{<:RaggedArray{Float32}},
    target_bundle::StrandPair{<:RaggedArray{Float32}},
    query_anchors::AnchorCSR,
    target_anchors::AnchorCSR,
    query_strand::Int,
    target_strand::Int,
    orientation_label::String,
    search_range::Int,
    window_radius::Int,
    realign_window::Int,
    metric::AbstractProfileMetric,
    min_logfpr::Float32=0.0f0,
)
    min_logfpr > 0.0f0 && return _score_orientation_pair_csr(
        query_bundle,
        target_bundle,
        query_anchors,
        target_anchors,
        query_strand,
        target_strand,
        orientation_label,
        search_range,
        window_radius,
        realign_window,
        metric,
    )

    query_scores = query_strand == 1 ? query_bundle.forward : query_bundle.reverse
    target_scores = target_strand == 1 ? target_bundle.forward : target_bundle.reverse
    best_score = 0.0f0
    best_shift = 0
    best_n_sites = 0
    for shift in (-search_range):search_range
        score, n_sites = _score_shift_best!(
            query_scores,
            target_scores,
            query_anchors,
            target_anchors,
            shift,
            window_radius,
            realign_window,
            metric,
        )
        if Float64(score) > Float64(best_score) || (
            Float64(score) == Float64(best_score) && (
                n_sites > best_n_sites ||
                (n_sites == best_n_sites && abs(shift) < abs(best_shift))
            )
        )
            best_score = score
            best_shift = shift
            best_n_sites = n_sites
        end
    end
    return (best_score, best_shift, best_n_sites, orientation_label)
end

"""
    ProfileConfig

Configuration for profile comparison.

Fields:
- `metric::AbstractProfileMetric`: profile metric type.
- `search_range::Int`: maximum shift to search (default 10).
- `window_radius::Int`: half-window size for site windows (default 10).
- `realign_window::Int`: realignment search radius (default 3).
- `min_logfpr::Float32`: minimum log FPR for threshold anchors (0 = best anchors).
"""
struct ProfileConfig{M<:AbstractProfileMetric}
    metric::M
    search_range::Int
    window_radius::Int
    realign_window::Int
    min_logfpr::Float32

    function ProfileConfig(
        metric::M,
        search_range::Int,
        window_radius::Int,
        realign_window::Int,
        min_logfpr::Float32,
    ) where {M<:AbstractProfileMetric}
        search_range >= 0 || throw(ArgumentError("search_range must be non-negative."))
        window_radius >= 0 || throw(ArgumentError("window_radius must be non-negative."))
        realign_window >= 0 || throw(ArgumentError("realign_window must be non-negative."))
        isfinite(min_logfpr) || throw(ArgumentError("min_logfpr must be finite."))
        return new{M}(metric, search_range, window_radius, realign_window, min_logfpr)
    end
end

function ProfileConfig(;
    metric::AbstractProfileMetric=OverlapCoefficient(),
    search_range::Int=10,
    window_radius::Int=10,
    realign_window::Int=3,
    min_logfpr::Real=0.0,
)
    return ProfileConfig(
        metric, search_range, window_radius, realign_window, Float32(min_logfpr)
    )
end

"""
    profile_compare(query_bundle, query_anchors, target_bundle, target_anchors, config::ProfileConfig)

Compare two normalized profile bundles using pre-collected anchor CSRs and
return
`(score::Float32, offset::Int, orientation::String, n_sites::Int, metric_name::String)`.

`query_anchors` and `target_anchors` are `(forward_csr, reverse_csr)` tuples.

Scores all four orientation pairs (`++`, `+-`, `-+`, `--`) and selects the best
with deterministic tie-breaking per ADR 0006: higher score wins, then
orientation priority `++ > +- > -+ > --`.
"""
function profile_compare(
    query_bundle::StrandPair{<:RaggedArray{Float32}},
    query_anchors::Tuple{AnchorCSR,AnchorCSR},
    target_bundle::StrandPair{<:RaggedArray{Float32}},
    target_anchors::Tuple{AnchorCSR,AnchorCSR},
    config::ProfileConfig,
)
    metric = config.metric

    # Score all four orientation pairs
    best_score = 0.0f0
    best_shift = 0
    best_n_sites = 0
    best_orientation = "++"
    best_rank = typemax(Int)

    for (label, q_strand, t_strand) in _orientation_pairs(query_bundle, target_bundle)
        qa = q_strand == 1 ? query_anchors[1] : query_anchors[2]
        ta = t_strand == 1 ? target_anchors[1] : target_anchors[2]
        result = _score_orientation_pair(
            query_bundle,
            target_bundle,
            qa,
            ta,
            q_strand,
            t_strand,
            label,
            config.search_range,
            config.window_radius,
            config.realign_window,
            metric,
            config.min_logfpr,
        )
        score, shift, n_sites, _ = result

        # Tie-breaking: score, site count, shift magnitude, then orientation rank.
        rank = PROFILE_ORIENTATION_RANK[label]
        if Float64(score) > Float64(best_score) || (
            Float64(score) == Float64(best_score) && (
                n_sites > best_n_sites || (
                    n_sites == best_n_sites && (
                        abs(shift) < abs(best_shift) ||
                        (abs(shift) == abs(best_shift) && rank < best_rank)
                    )
                )
            )
        )
            best_score = score
            best_shift = shift
            best_n_sites = n_sites
            best_orientation = label
            best_rank = rank
        end
    end

    return (best_score, best_shift, best_orientation, best_n_sites, metric_name(metric))
end

"""
    profile_compare(query_bundle, target_bundle, config::ProfileConfig)

Compare two normalized profile bundles and return
`(score::Float32, offset::Int, orientation::String, n_sites::Int, metric_name::String)`.

Collects anchors for both bundles internally. For one-to-many comparison where
the query anchors should be reused, use the variant with pre-computed anchors.

Scores all four orientation pairs (`++`, `+-`, `-+`, `--`) and selects the best
with deterministic tie-breaking per ADR 0006: higher score wins, then
orientation priority `++ > +- > -+ > --`.
"""
function profile_compare(
    query_bundle::StrandPair{<:RaggedArray{Float32}},
    target_bundle::StrandPair{<:RaggedArray{Float32}},
    config::ProfileConfig,
)
    threshold = config.min_logfpr
    query_anchors = _collect_both_anchors(query_bundle, threshold)
    target_anchors = _collect_both_anchors(target_bundle, threshold)
    return profile_compare(
        query_bundle, query_anchors, target_bundle, target_anchors, config
    )
end

# ── Prepared profile (one-to-many reuse) ───────────────────────────────────────

"""
    PreparedProfile

A pre-normalized profile bundle with pre-collected anchors, ready for
repeated comparison against multiple targets without re-computing the
query side.

Fields:
- `name::String`: profile name.
- `bundle::StrandPair{T}`: normalized scores.
- `anchors::Tuple{AnchorCSR,AnchorCSR}`: `(forward, reverse)` anchor CSRs.
"""
struct PreparedProfile{T}
    name::String
    bundle::T
    anchors::Tuple{AnchorCSR,AnchorCSR}
    min_logfpr::Float32
end

"""
    modelname(profile::PreparedProfile)

Return the prepared profile name. `PreparedProfile` is an
`AbstractProfileSource` and uses `modelname` as its public name accessor.
"""
modelname(profile::PreparedProfile) = profile.name

function PreparedProfile(
    name::String,
    bundle::StrandPair{<:RaggedArray{Float32}},
    anchors::Tuple{AnchorCSR,AnchorCSR},
    min_logfpr::Float32,
)
    threshold = min_logfpr
    isfinite(threshold) || throw(ArgumentError("min_logfpr must be finite."))
    n_rows = nrows(bundle.forward)
    n_rows == nrows(bundle.reverse) ||
        throw(ArgumentError("prepared strand bundles must have equal row counts."))
    length(anchors[1].offsets) == n_rows + 1 ||
        throw(ArgumentError("forward anchor rows do not match the profile bundle."))
    length(anchors[2].offsets) == n_rows + 1 ||
        throw(ArgumentError("reverse anchor rows do not match the profile bundle."))
    for (csr, strand) in zip(anchors, (bundle.forward, bundle.reverse))
        for row_index in 1:n_rows
            for position in
                csr.positions[csr.offsets[row_index]:(csr.offsets[row_index + 1] - 1)]
                1 <= position <= rowlength(strand, row_index) ||
                    throw(ArgumentError("anchor position is outside its profile row."))
            end
        end
    end
    return PreparedProfile{typeof(bundle)}(name, bundle, anchors, threshold)
end

function PreparedProfile(
    name::String,
    bundle::StrandPair{<:RaggedArray{Float32}},
    anchors::Tuple{AnchorCSR,AnchorCSR},
    min_logfpr::Real,
)
    return PreparedProfile(name, bundle, anchors, Float32(min_logfpr))
end

function PreparedProfile(
    name::String,
    bundle::StrandPair{<:RaggedArray{Float32}},
    anchors::Tuple{AnchorCSR,AnchorCSR},
)
    return PreparedProfile(name, bundle, anchors, 0.0f0)
end

function PreparedProfile(
    name::AbstractString,
    bundle::StrandPair{<:RaggedArray{Float32}},
    anchors::Tuple{AnchorCSR,AnchorCSR},
    min_logfpr::Real,
)
    return PreparedProfile(String(name), bundle, anchors, min_logfpr)
end

function PreparedProfile(
    name::AbstractString,
    bundle::StrandPair{<:RaggedArray{Float32}},
    anchors::Tuple{AnchorCSR,AnchorCSR},
)
    return PreparedProfile(String(name), bundle, anchors)
end

"""
    _collect_both_anchors(bundle::StrandPair{<:RaggedArray{Float32}}, threshold::Float32)

Collect anchors for both strands and return `(forward_csr, reverse_csr)`.
"""
function _collect_both_anchors(
    bundle::StrandPair{<:RaggedArray{Float32}}, threshold::Float32
)
    n_rows = nrows(bundle.forward)
    fwd_rows, fwd_pos = collect_anchors(bundle.forward, threshold)
    fwd_csr = build_anchor_csr(fwd_rows, fwd_pos, n_rows)
    if bundle.forward === bundle.reverse
        return (fwd_csr, fwd_csr)
    end
    rev_rows, rev_pos = collect_anchors(bundle.reverse, threshold)
    return (fwd_csr, build_anchor_csr(rev_rows, rev_pos, n_rows))
end

"""
    prepare_profile(model::ScoreProfile; min_logfpr::Float32=0.0f0)

Prepare a [`ScoreProfile`](@ref) for repeated comparison: fit normalization
from the profile's own scores, apply it, and collect anchors for both
strands. Returns a [`PreparedProfile`](@ref).

Keyword arguments:
- `min_logfpr::Float32=0.0`: minimum log FPR for threshold anchors (0 = best anchors).
"""
function prepare_profile(
    model::ScoreProfile; min_logfpr::Real=0.0, execution::ExecutionPolicy=SerialExecution()
)
    threshold = Float32(min_logfpr)
    _, normalized = _fit_transform_empirical(model.scores)
    norm_bundle = StrandPair(normalized, normalized)
    anchors = _collect_both_anchors(norm_bundle, threshold)
    return PreparedProfile(String(modelname(model)), norm_bundle, anchors, threshold)
end

"""
    prepare_profile(model::AbstractMotifModel, sequences::EncodedSequenceBatch;
                    background=nothing, min_logfpr=0.0f0)

Scan a motif model against `sequences` to produce a normalized profile bundle
with pre-collected anchors, ready for repeated comparison. The normalization
is fitted from `background` (falls back to `sequences` when `background=nothing`).

Returns a [`PreparedProfile`](@ref).
"""
function prepare_profile(
    model::AbstractMotifModel,
    sequences::EncodedSequenceBatch;
    background::Union{EncodedSequenceBatch,Nothing}=nothing,
    min_logfpr::Real=0.0,
    execution::ExecutionPolicy=SerialExecution(),
)
    threshold = Float32(min_logfpr)
    validate_model(model; capability=:compare)
    raw = _scan_model_batch(model, sequences; strands=BothStrands(), execution=execution)
    bg = background === nothing ? sequences : background
    if bg === sequences
        _, norm_bundle = _fit_transform_empirical(raw)
    else
        bg_raw = _scan_model_batch(model, bg; strands=BothStrands(), execution=execution)
        table = fit(EmpiricalLogTail(), flatten_bundle(bg_raw))
        norm_bundle = normalize_bundle(table, raw)
    end
    anchors = _collect_both_anchors(norm_bundle, threshold)
    return PreparedProfile(String(modelname(model)), norm_bundle, anchors, threshold)
end

# ── PreparedProfile compare methods ────────────────────────────────────────────

"""
    compare(query::PreparedProfile, target::ScoreProfile; metric=:co, kwargs...)

Compare a [`PreparedProfile`](@ref) (pre-normalized, pre-anchored) against a
[`ScoreProfile`](@ref) target. The query normalization and anchor collection
are reused; only the target is prepared.

Keyword arguments:
- `metric`: profile metric (`:co`, `:co_rowwise`, `:dice`, `:dice_rowwise`,
  `:cosine`, or a typed `AbstractProfileMetric`). Default `:co`.
- `search_range::Int=10`, `window_radius::Int=10`, `realign_window::Int=3`,
  `min_logfpr::Float32=0.0`.
"""
function compare(
    query::PreparedProfile,
    target::ScoreProfile;
    metric::Union{AbstractString,Symbol,AbstractProfileMetric}=:co,
    search_range::Int=10,
    window_radius::Int=10,
    realign_window::Int=3,
    min_logfpr::Union{Nothing,Real}=nothing,
)
    m = _resolve_profile_metric(metric)
    threshold = min_logfpr === nothing ? query.min_logfpr : Float32(min_logfpr)
    threshold == query.min_logfpr ||
        throw(ArgumentError("min_logfpr differs from the prepared query threshold."))
    _, target_scores = _fit_transform_empirical(target.scores)
    target_norm = StrandPair(target_scores, target_scores)
    target_anchors = _collect_both_anchors(target_norm, threshold)
    config = ProfileConfig(;
        metric=m,
        search_range=search_range,
        window_radius=window_radius,
        realign_window=realign_window,
        min_logfpr=threshold,
    )
    score, shift, orientation, n_sites, metric_str = profile_compare(
        query.bundle, query.anchors, target_norm, target_anchors, config
    )
    return ComparisonResult(
        modelname(query), modelname(target), score, shift, orientation, metric_str, n_sites
    )
end

function compare(
    query::PreparedProfile,
    target::PreparedProfile;
    metric::Union{AbstractString,Symbol,AbstractProfileMetric}=:co,
    search_range::Int=10,
    window_radius::Int=10,
    realign_window::Int=3,
)
    query.min_logfpr == target.min_logfpr ||
        throw(ArgumentError("prepared profiles use different min_logfpr thresholds."))
    config = ProfileConfig(;
        metric=_resolve_profile_metric(metric),
        search_range,
        window_radius,
        realign_window,
        min_logfpr=query.min_logfpr,
    )
    score, shift, orientation, n_sites, metric_str = profile_compare(
        query.bundle, query.anchors, target.bundle, target.anchors, config
    )
    return ComparisonResult(
        modelname(query), modelname(target), score, shift, orientation, metric_str, n_sites
    )
end

"""
    compare(query::PreparedProfile, targets::Vector{ScoreProfile}; kwargs...)

One-to-many comparison: compare a single prepared query against multiple
targets, reusing the query's normalized bundle and pre-collected anchors.
Returns a `Vector{ComparisonResult}`.

Each target is prepared (normalized, anchors collected) independently.
"""
function compare(
    query::PreparedProfile,
    targets::AbstractVector{<:ScoreProfile};
    execution::ExecutionPolicy=SerialExecution(),
    metric::Union{AbstractString,Symbol,AbstractProfileMetric}=:co,
    search_range::Int=10,
    window_radius::Int=10,
    realign_window::Int=3,
    min_logfpr::Union{Nothing,Real}=nothing,
)
    threshold = min_logfpr === nothing ? query.min_logfpr : Float32(min_logfpr)
    threshold == query.min_logfpr ||
        throw(ArgumentError("min_logfpr differs from the prepared query threshold."))
    m = _resolve_profile_metric(metric)
    config = ProfileConfig(;
        metric=m, search_range, window_radius, realign_window, min_logfpr=threshold
    )
    results = Vector{ComparisonResult}(undef, length(targets))
    isempty(targets) && return results
    _parallel_for(execution, length(targets)) do i
        target = targets[i]
        _, target_scores = _fit_transform_empirical(target.scores)
        target_bundle = StrandPair(target_scores, target_scores)
        target_anchors = _collect_both_anchors(target_bundle, threshold)
        score, shift, orientation, n_sites, metric_str = profile_compare(
            query.bundle, query.anchors, target_bundle, target_anchors, config
        )
        return results[i] = ComparisonResult(
            modelname(query),
            modelname(target),
            score,
            shift,
            orientation,
            metric_str,
            n_sites,
        )
    end
    return results
end

function compare(
    query::PreparedProfile,
    targets::AbstractVector{<:PreparedProfile};
    execution::ExecutionPolicy=SerialExecution(),
    metric::Union{AbstractString,Symbol,AbstractProfileMetric}=:co,
    search_range::Int=10,
    window_radius::Int=10,
    realign_window::Int=3,
)
    config = ProfileConfig(;
        metric=_resolve_profile_metric(metric),
        search_range,
        window_radius,
        realign_window,
        min_logfpr=query.min_logfpr,
    )
    results = Vector{ComparisonResult}(undef, length(targets))
    isempty(targets) && return results
    _parallel_for(execution, length(targets)) do i
        target = targets[i]
        query.min_logfpr == target.min_logfpr ||
            throw(ArgumentError("prepared profiles use different min_logfpr thresholds."))
        score, shift, orientation, n_sites, metric_str = profile_compare(
            query.bundle, query.anchors, target.bundle, target.anchors, config
        )
        return results[i] = ComparisonResult(
            modelname(query),
            modelname(target),
            score,
            shift,
            orientation,
            metric_str,
            n_sites,
        )
    end
    return results
end

"""
    compare(query::PreparedProfile, target::AbstractMotifModel,
            sequences::EncodedSequenceBatch; metric=:co,
            execution=SerialExecution(), kwargs...)

Compare a [`PreparedProfile`](@ref) against a motif model target by scanning
the target against `sequences` and comparing profiles. The query's normalized
bundle and anchors are reused.
"""
function compare(
    query::PreparedProfile,
    target::AbstractMotifModel,
    sequences::EncodedSequenceBatch;
    metric::Union{AbstractString,Symbol,AbstractProfileMetric}=:co,
    search_range::Int=10,
    window_radius::Int=10,
    realign_window::Int=3,
    min_logfpr::Union{Nothing,Real}=nothing,
    background::Union{EncodedSequenceBatch,Nothing}=nothing,
    execution::ExecutionPolicy=SerialExecution(),
)
    m = _resolve_profile_metric(metric)
    threshold = min_logfpr === nothing ? query.min_logfpr : Float32(min_logfpr)
    threshold == query.min_logfpr ||
        throw(ArgumentError("min_logfpr differs from the prepared query threshold."))
    target_norm = _resolve_profile_bundle(
        target, sequences, background; execution=execution
    )
    target_anchors = _collect_both_anchors(target_norm, threshold)
    config = ProfileConfig(;
        metric=m,
        search_range=search_range,
        window_radius=window_radius,
        realign_window=realign_window,
        min_logfpr=threshold,
    )
    score, shift, orientation, n_sites, metric_str = profile_compare(
        query.bundle, query.anchors, target_norm, target_anchors, config
    )
    return ComparisonResult(
        modelname(query), modelname(target), score, shift, orientation, metric_str, n_sites
    )
end

"""
    compare(query::AbstractMotifModel, target::PreparedProfile,
            sequences::EncodedSequenceBatch; metric=:co,
            execution=SerialExecution(), kwargs...)

Compare a motif model query (scanned against `sequences`) against a
[`PreparedProfile`](@ref) target. The target's normalized bundle and anchors
are reused.
"""
function compare(
    query::AbstractMotifModel,
    target::PreparedProfile,
    sequences::EncodedSequenceBatch;
    metric::Union{AbstractString,Symbol,AbstractProfileMetric}=:co,
    search_range::Int=10,
    window_radius::Int=10,
    realign_window::Int=3,
    min_logfpr::Union{Nothing,Real}=nothing,
    background::Union{EncodedSequenceBatch,Nothing}=nothing,
    execution::ExecutionPolicy=SerialExecution(),
)
    m = _resolve_profile_metric(metric)
    threshold = min_logfpr === nothing ? target.min_logfpr : Float32(min_logfpr)
    threshold == target.min_logfpr ||
        throw(ArgumentError("min_logfpr differs from the prepared target threshold."))
    query_norm = _resolve_profile_bundle(query, sequences, background; execution=execution)
    query_anchors = _collect_both_anchors(query_norm, threshold)
    config = ProfileConfig(;
        metric=m,
        search_range=search_range,
        window_radius=window_radius,
        realign_window=realign_window,
        min_logfpr=threshold,
    )
    score, shift, orientation, n_sites, metric_str = profile_compare(
        query_norm, query_anchors, target.bundle, target.anchors, config
    )
    return ComparisonResult(
        modelname(query), modelname(target), score, shift, orientation, metric_str, n_sites
    )
end

"""
    compare(query::PreparedProfile, targets, sequences; execution=...)

Compare one prepared query against motif-model targets. Each target owns a
serial scan, normalization, anchor collection, and alignment path.
"""
function compare(
    query::PreparedProfile,
    targets::AbstractVector{<:AbstractMotifModel},
    sequences::EncodedSequenceBatch;
    execution::ExecutionPolicy=SerialExecution(),
    metric::Union{AbstractString,Symbol,AbstractProfileMetric}=:co,
    search_range::Int=10,
    window_radius::Int=10,
    realign_window::Int=3,
    min_logfpr::Union{Nothing,Real}=nothing,
    background::Union{EncodedSequenceBatch,Nothing}=nothing,
)
    threshold = min_logfpr === nothing ? query.min_logfpr : Float32(min_logfpr)
    threshold == query.min_logfpr ||
        throw(ArgumentError("min_logfpr differs from the prepared query threshold."))
    config = ProfileConfig(;
        metric=_resolve_profile_metric(metric),
        search_range,
        window_radius,
        realign_window,
        min_logfpr=threshold,
    )
    results = Vector{ComparisonResult}(undef, length(targets))
    isempty(targets) && return results
    _parallel_for(execution, length(targets)) do i
        target = targets[i]
        target_bundle = _resolve_profile_bundle(
            target, sequences, background; execution=SerialExecution()
        )
        target_anchors = _collect_both_anchors(target_bundle, threshold)
        score, shift, orientation, n_sites, metric_str = profile_compare(
            query.bundle, query.anchors, target_bundle, target_anchors, config
        )
        return results[i] = ComparisonResult(
            modelname(query),
            modelname(target),
            score,
            shift,
            orientation,
            metric_str,
            n_sites,
        )
    end
    return results
end

"""
    compare(query_model, targets, sequences; execution=...)

Prepare the query once, then compare it against motif-model targets at the
outer target level. Inner target work is explicitly serial.
"""
function compare(
    query::AbstractMotifModel,
    targets::AbstractVector{<:AbstractMotifModel},
    sequences::EncodedSequenceBatch;
    execution::ExecutionPolicy=SerialExecution(),
    metric::Union{AbstractString,Symbol,AbstractProfileMetric}=:co,
    search_range::Int=10,
    window_radius::Int=10,
    realign_window::Int=3,
    min_logfpr::Real=0.0,
    background::Union{EncodedSequenceBatch,Nothing}=nothing,
)
    prepared_query = prepare_profile(query, sequences; background, min_logfpr, execution)
    return compare(
        prepared_query,
        targets,
        sequences;
        execution,
        metric,
        search_range,
        window_radius,
        realign_window,
        background,
    )
end
