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
2. Fit hybrid empirical log-tail normalization from each model's own scores.
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
    normalization::AbstractNormalizationStrategy=HybridEmpiricalLogTail(),
    execution::Execution=Execution(),
    cache=nothing,
)
    m = _resolve_profile_metric(metric)

    return compare(
        prepare_profile(
            query; min_logfpr=Float32(min_logfpr), normalization, execution, cache
        ),
        prepare_profile(
            target; min_logfpr=Float32(min_logfpr), normalization, execution, cache
        );
        metric=m,
        search_range=search_range,
        window_radius=window_radius,
        realign_window=realign_window,
        execution,
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
    normalization::AbstractNormalizationStrategy=HybridEmpiricalLogTail(),
    execution::Execution=Execution(),
    cache=nothing,
)
    m = _resolve_profile_metric(metric)
    threshold = Float32(min_logfpr)
    prepared_query = prepare_profile(
        query,
        sequences;
        background=background,
        min_logfpr=threshold,
        normalization=normalization,
        execution=execution,
        cache=cache,
    )
    prepared_target = prepare_profile(
        target,
        sequences;
        background=background,
        min_logfpr=threshold,
        normalization=normalization,
        execution=execution,
        cache=cache,
    )
    config = ProfileConfig(;
        metric=m,
        search_range=search_range,
        window_radius=window_radius,
        realign_window=realign_window,
        min_logfpr=threshold,
    )
    score, shift, orientation, n_sites, metric_str = profile_compare(
        prepared_query.bundle,
        prepared_query.anchors,
        prepared_target.bundle,
        prepared_target.anchors,
        config;
        execution,
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
