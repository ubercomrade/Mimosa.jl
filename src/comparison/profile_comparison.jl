# Profile comparison dispatch: compare ScoreProfile models using the profile algorithm.

"""
    compare(query::ScoreProfile, target::ScoreProfile; metric=:co, kwargs...)

Compare two [`ScoreProfile`](@ref) models using the window-based profile
comparison algorithm. Returns a [`ComparisonResult`](@ref) with deterministic
tie-breaking per ADR 0006.

Keyword arguments:
- `metric`: profile metric (`:co`, `:co_rowwise`, `:dice`, `:dice_rowwise`,
  `:cosine`, or a typed `AbstractProfileMetric`). Default `:co`.
- `search_range::Int=10`: maximum shift to search.
- `window_radius::Int=10`: half-window size for site windows.
- `realign_window::Int=3`: realignment search radius for target anchors.
- `min_logfpr::Float32=0.0`: minimum log FPR for threshold anchors (0 = best anchors).

The comparison pipeline:
1. Resolve profile bundles (both strands = same scores for ScoreProfile).
2. Fit `EmpiricalLogTail` normalization from each model's own scores.
3. Apply normalization to both strands.
4. Collect anchors (best per row or threshold).
5. Score all four orientation pairs across all shifts.
6. Select best with deterministic tie-breaking.
"""
function compare(
    query::ScoreProfile,
    target::ScoreProfile;
    metric::Union{AbstractString,Symbol,AbstractProfileMetric}=:co,
    search_range::Int=10,
    window_radius::Int=10,
    realign_window::Int=3,
    min_logfpr::Real=0.0,
)
    m = _resolve_profile_metric(metric)

    return compare(
        prepare_profile(query; min_logfpr=Float32(min_logfpr)),
        prepare_profile(target; min_logfpr=Float32(min_logfpr));
        metric=m,
        search_range=search_range,
        window_radius=window_radius,
        realign_window=realign_window,
    )
end

function compare(
    query::ScoreProfile,
    target::AbstractMotifModel,
    sequences::EncodedSequenceBatch;
    kwargs...,
)
    return throw(
        ArgumentError(
            "mixed ScoreProfile/motif comparison is unsupported; prepare both inputs as profiles first.",
        ),
    )
end

function compare(
    query::AbstractMotifModel,
    target::ScoreProfile,
    sequences::EncodedSequenceBatch;
    kwargs...,
)
    return throw(
        ArgumentError(
            "mixed motif/ScoreProfile comparison is unsupported; prepare both inputs as profiles first.",
        ),
    )
end

function compare(
    query::ScoreProfile, target::ScoreProfile, sequences::EncodedSequenceBatch; kwargs...
)
    return throw(ArgumentError("ScoreProfile comparison does not consume sequences."))
end

# ── Motif-derived profile comparison ──────────────────────────────────────────

"""
    _resolve_profile_bundle(model, sequences, background_sequences; kwargs...)

Scan a motif model against `sequences` to produce a raw strand-aware profile
bundle, then fit and apply `EmpiricalLogTail` normalization from
`background_sequences` (falls back to `sequences` when `background_sequences`
is `nothing`).

Returns `StrandPair{RaggedArray{Float32}}` with normalized scores.
"""
function _resolve_profile_bundle(
    model::AbstractMotifModel,
    sequences::EncodedSequenceBatch,
    background_sequences::Union{EncodedSequenceBatch,Nothing};
    execution::ExecutionPolicy=SerialExecution(),
)
    validate_model(model; capability=:compare)
    raw = _scan_model_batch(model, sequences; strands=BothStrands(), execution=execution)
    bg = background_sequences === nothing ? sequences : background_sequences
    if bg === sequences
        _, normalized = _fit_transform_empirical(raw)
        return normalized
    end
    bg_raw = _scan_model_batch(model, bg; strands=BothStrands(), execution=execution)
    table = fit(EmpiricalLogTail(), flatten_bundle(bg_raw))
    return normalize_bundle(table, raw)
end

"""
    compare(query::AbstractMotifModel, target::AbstractMotifModel,
            sequences::EncodedSequenceBatch; metric=:co, kwargs...)

Compare two motif models via the profile-based comparison strategy: scan both
models against `sequences` to produce score profiles, normalize, and compare
using the window-based profile algorithm.

This is the Julia equivalent of Python's `strategy_profile`.

Keyword arguments:
- `metric`: profile metric (`:co`, `:co_rowwise`, `:dice`, `:dice_rowwise`,
  `:cosine`, or a typed `AbstractProfileMetric`). Default `:co`.
- `search_range::Int=10`, `window_radius::Int=10`, `realign_window::Int=3`,
  `min_logfpr::Float32=0.0`.
- `background::Union{EncodedSequenceBatch,Nothing}=nothing`: optional separate
  background sequences for normalization. Falls back to `sequences`.
"""
function compare(
    query::AbstractMotifModel,
    target::AbstractMotifModel,
    sequences::EncodedSequenceBatch;
    metric::Union{AbstractString,Symbol,AbstractProfileMetric}=:co,
    search_range::Int=10,
    window_radius::Int=10,
    realign_window::Int=3,
    min_logfpr::Real=0.0,
    background::Union{EncodedSequenceBatch,Nothing}=nothing,
    execution::ExecutionPolicy=SerialExecution(),
)
    m = _resolve_profile_metric(metric)
    threshold = Float32(min_logfpr)
    query_norm = _resolve_profile_bundle(query, sequences, background; execution=execution)
    target_norm = _resolve_profile_bundle(
        target, sequences, background; execution=execution
    )
    config = ProfileConfig(;
        metric=m,
        search_range=search_range,
        window_radius=window_radius,
        realign_window=realign_window,
        min_logfpr=threshold,
    )
    query_anchors = _collect_both_anchors(query_norm, threshold)
    target_anchors = _collect_both_anchors(target_norm, threshold)
    score, shift, orientation, n_sites, metric_str = profile_compare(
        query_norm, query_anchors, target_norm, target_anchors, config
    )
    return ComparisonResult(
        String(modelname(query)),
        String(modelname(target)),
        score,
        shift,
        orientation,
        metric_str,
        n_sites,
    )
end
