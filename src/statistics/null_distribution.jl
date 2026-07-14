# Null distribution type, build workflow, and result annotation.

"""
    NullPair

A single contributing comparison pair in a null distribution.

Fields:
- `query::String`: query model name.
- `target::String`: target model name.
- `score::Float64`: raw comparison score.
"""
struct NullPair
    query::String
    target::String
    score::Float64
end

"""Canonical numerical and data provenance required to annotate a null."""
struct ProfileComparisonContract
    metric::String
    search_range::Int
    window_radius::Int
    realign_window::Int
    min_logfpr::Float32
    normalization_version::String
    alignment_version::String
    sequence_fingerprint::String
    background_fingerprint::String
    raw_scores_fingerprint::String
end

"""
    NullDistribution

A fitted null distribution for significance testing.

Fields:
- `strategy::String`: fixed profile marker (`"profile"`).
- `metric::String`: comparison metric name.
- `fit::GEVFitResult`: fitted GEV distribution or failure.
- `raw_scores::Vector{Float64}`: pooled raw comparison scores.
- `pairs::Vector{NullPair}`: contributing comparison pairs.
- `n_null::Int`: number of null comparisons.
- `n_queries::Int`: number of queries used.
- `skipped::Vector{NamedTuple{(:query, :reason),Tuple{String,String}}}`: skipped queries.
- `model_collection_fingerprint::Union{Nothing,String}`: SHA-256 of model collection.
- `relation_fingerprint::Union{Nothing,String}`: SHA-256 of relation data.
- `sequence_fingerprint::String`: fingerprint of sequences (or `"none"`).
- `background_fingerprint::String`: fingerprint of background (or `"none"`).
"""
struct NullDistribution
    strategy::String
    metric::String
    fit::GEVFitResult
    raw_scores::Vector{Float64}
    pairs::Vector{NullPair}
    n_null::Int
    n_queries::Int
    skipped::Vector{NamedTuple{(:query, :reason),Tuple{String,String}}}
    model_collection_fingerprint::Union{Nothing,String}
    relation_fingerprint::Union{Nothing,String}
    sequence_fingerprint::String
    background_fingerprint::String
    contract::ProfileComparisonContract
end

function NullDistribution(
    strategy::String,
    metric::String,
    fit::GEVFitResult,
    raw_scores::Vector{Float64},
    pairs::Vector{NullPair},
    n_null::Int,
    n_queries::Int,
    skipped,
    model_collection_fingerprint,
    relation_fingerprint,
    sequence_fingerprint::String,
    background_fingerprint::String,
)
    contract = ProfileComparisonContract(
        metric,
        10,
        10,
        3,
        0.0f0,
        "empirical-log-tail-v1",
        "profile-alignment-v1",
        sequence_fingerprint,
        background_fingerprint,
        content_fingerprint(raw_scores),
    )
    return NullDistribution(
        strategy,
        metric,
        fit,
        raw_scores,
        pairs,
        n_null,
        n_queries,
        skipped,
        model_collection_fingerprint,
        relation_fingerprint,
        sequence_fingerprint,
        background_fingerprint,
        contract,
    )
end

"""
    NullBuildConfig

Configuration for building a null distribution.

Fields:
- `metric`: comparison metric (Symbol, String, or typed metric).
- `min_null_targets::Int`: minimum eligible targets required per query.
- `strict::Bool`: if `true`, raise an error when a query has too few targets.
"""
struct NullBuildConfig{M<:AbstractProfileMetric}
    metric::M
    min_null_targets::Int
    strict::Bool

    function NullBuildConfig(
        metric::M, min_null_targets::Int, strict::Bool
    ) where {M<:AbstractProfileMetric}
        min_null_targets > 0 || throw(ArgumentError("min_null_targets must be positive."))
        return new{M}(metric, min_null_targets, strict)
    end
end

function NullBuildConfig(; metric=nothing, min_null_targets::Int=1, strict::Bool=false)
    resolved_metric = _resolve_profile_metric(isnothing(metric) ? :co : metric)
    return NullBuildConfig(resolved_metric, min_null_targets, strict)
end

"""
    NullBuildResult

Result of building a null distribution, including the distribution and build
statistics.
"""
struct NullBuildResult
    distribution::NullDistribution
    total_comparisons::Int
end

"""
    build_null(models, relations; sequences=batch, metric=:co, min_null_targets=1,
               strict=false, execution=SerialExecution(), kwargs...)

Build a pooled null distribution from all eligible query-target comparisons.

For each model (query), compares it against all eligible targets (motifs from
a different group) and collects raw scores. The pooled scores are fitted to a
GEV distribution.

# Arguments
- `models::AbstractVector`: vector of motif models (e.g. `PWM`).
- `relations::GroupRelations`: group relations from [`parse_group_relations`](@ref).
- `metric`: comparison metric.
- `min_null_targets`: minimum eligible targets per query (default 1).
- `strict`: if `true`, raise an error when a query has too few targets.
- `execution`: [`ExecutionPolicy`](@ref) for parallel comparison of query-target
  pairs. Default `SerialExecution()`.
- `kwargs...`: additional keyword arguments passed to `compare`.

Returns a [`NullBuildResult`](@ref).

Under `ThreadedExecution`, comparisons are processed in parallel at the
top level. Results are collected into pre-allocated slots indexed by
the original comparison order, so the pooled score order and fit are
identical to `SerialExecution`.
"""
function build_null(
    models::AbstractVector,
    relations::GroupRelations;
    metric=nothing,
    min_null_targets::Int=1,
    strict::Bool=false,
    execution::ExecutionPolicy=SerialExecution(),
    sequences::EncodedSequenceBatch,
    background::Union{Nothing,EncodedSequenceBatch}=nothing,
    search_range::Int=10,
    window_radius::Int=10,
    realign_window::Int=3,
    min_logfpr::Real=0.0,
    kwargs...,
)
    config = NullBuildConfig(;
        metric=metric, min_null_targets=min_null_targets, strict=strict
    )
    return build_null(
        models,
        relations,
        config;
        execution=execution,
        sequences=sequences,
        background=background,
        search_range=search_range,
        window_radius=window_radius,
        realign_window=realign_window,
        min_logfpr=min_logfpr,
        kwargs...,
    )
end

function build_null(
    models::AbstractVector,
    relations::GroupRelations,
    config::NullBuildConfig{<:AbstractProfileMetric};
    execution::ExecutionPolicy=SerialExecution(),
    sequences::EncodedSequenceBatch,
    background::Union{Nothing,EncodedSequenceBatch}=nothing,
    search_range::Int=10,
    window_radius::Int=10,
    realign_window::Int=3,
    min_logfpr::Real=0.0,
    kwargs...,
)
    isempty(kwargs) || throw(ArgumentError("unsupported profile null-build option."))
    for model in models
        model isa AbstractMotifModel && validate_model(model; capability=:cache)
    end
    names = String[modelname(model) for model in models]
    length(unique(names)) == length(names) ||
        throw(ArgumentError("model names must be unique for null construction."))
    all(!isempty, names) || throw(ArgumentError("model names must not be empty."))
    search_range >= 0 || throw(ArgumentError("search_range must be non-negative."))
    window_radius >= 0 || throw(ArgumentError("window_radius must be non-negative."))
    realign_window >= 0 || throw(ArgumentError("realign_window must be non-negative."))
    isfinite(min_logfpr) || throw(ArgumentError("min_logfpr must be finite."))
    prepared = Vector{PreparedProfile}(undef, length(models))
    _parallel_for(execution, length(models)) do i
        model = models[i]
        return prepared[i] = if model isa ScoreProfile
            prepare_profile(model; min_logfpr=min_logfpr)
        else
            prepare_profile(
                model,
                sequences;
                background=background,
                min_logfpr=min_logfpr,
                execution=SerialExecution(),
            )
        end
    end
    prepared_by_name = Dict(modelname(p) => p for p in prepared)
    result = _build_null(models, relations, config, execution) do q, t
        return compare(
            prepared_by_name[modelname(q)],
            prepared_by_name[modelname(t)];
            metric=config.metric,
            search_range=search_range,
            window_radius=window_radius,
            realign_window=realign_window,
        )
    end
    dist = result.distribution
    profiled_dist = NullDistribution(
        dist.strategy,
        dist.metric,
        dist.fit,
        dist.raw_scores,
        dist.pairs,
        dist.n_null,
        dist.n_queries,
        dist.skipped,
        dist.model_collection_fingerprint,
        dist.relation_fingerprint,
        sequence_fingerprint(sequences),
        isnothing(background) ? "none" : sequence_fingerprint(background),
        ProfileComparisonContract(
            dist.metric,
            search_range,
            window_radius,
            realign_window,
            Float32(min_logfpr),
            "empirical-log-tail-v1",
            "profile-alignment-v1",
            sequence_fingerprint(sequences),
            isnothing(background) ? "none" : sequence_fingerprint(background),
            content_fingerprint(dist.raw_scores),
        ),
    )
    return NullBuildResult(profiled_dist, result.total_comparisons)
end

function _build_null(
    compare_pair,
    models::AbstractVector,
    relations::GroupRelations,
    config::NullBuildConfig,
    execution::ExecutionPolicy,
)
    all(source -> source isa AbstractProfileSource, models) ||
        throw(ArgumentError("models must contain only AbstractProfileSource values."))
    by_name = Dict{String,AbstractProfileSource}()
    for model in models
        by_name[modelname(model)] = model
    end

    # Build the work schedule: list of (query, target) pairs to compare.
    #
    # Type note: work_pairs uses AbstractProfileSource values because the
    # collection may mix motif models and precomputed profiles.
    # This abstract element type is unavoidable without requiring homogeneous
    # model collections or a separate build_null method per concrete model type.
    # The inner comparison loop dispatches through the compare_pair closure,
    # which is type-stable per individual call. This is NOT a hot path — the
    # expensive work is inside compare_pair, not in the work_pairs iteration.
    # See PLAN_2.md E3 for the type/allocation audit rationale.
    work_pairs = Tuple{AbstractProfileSource,AbstractProfileSource}[]
    skipped = NamedTuple{(:query, :reason),Tuple{String,String}}[]
    n_queries = 0

    for query in models
        qname = modelname(query)
        target_names = eligible_targets(relations, qname)
        # Filter to known models and exclude self
        target_names = filter(n -> n != qname && haskey(by_name, n), target_names)

        if length(target_names) < config.min_null_targets
            reason = "only $(length(target_names)) null target(s); required $(config.min_null_targets)"
            push!(skipped, (query=qname, reason=reason))
            if config.strict
                throw(ArgumentError("Skipping null contribution for $qname: $reason"))
            end
            continue
        end

        n_queries += 1
        for target_name in target_names
            push!(work_pairs, (query, by_name[target_name]))
        end
    end

    total_comparisons = length(work_pairs)
    if total_comparisons == 0
        throw(
            ArgumentError(
                "Cannot build a null distribution: no eligible query-target comparisons were found.",
            ),
        )
    end

    # Pre-allocate result slots for raw scores
    raw_scores = Vector{Float64}(undef, total_comparisons)
    pairs = Vector{NullPair}(undef, total_comparisons)

    # Compare all pairs (serial or threaded)
    _parallel_for(execution, total_comparisons) do i
        q, t = work_pairs[i]
        result = compare_pair(q, t)
        score = Float64(result.score)
        raw_scores[i] = score
        return pairs[i] = NullPair(String(modelname(q)), String(modelname(t)), score)
    end

    fit_result = fit_gev(raw_scores)

    dist = NullDistribution(
        "profile",
        metric_name(config.metric),
        fit_result,
        raw_scores,
        pairs,
        length(raw_scores),
        n_queries,
        skipped,
        model_collection_fingerprint(AbstractProfileSource[models...]),
        _relation_fingerprint(relations),
        "none",
        "none",
    )

    return NullBuildResult(dist, total_comparisons)
end

function _relation_fingerprint(relations::GroupRelations)
    entries = String[]
    for name in sort!(collect(keys(relations.groups)))
        push!(entries, "$(name)=$(relations.groups[name])")
    end
    return content_fingerprint(join(entries, "\n"))
end

# ---------------------------------------------------------------------------
# Result annotation
# ---------------------------------------------------------------------------

"""
    AnnotatedResult

A comparison result enriched with significance values from a null distribution.

Fields (same as `ComparisonResult` plus significance):
- `query::String`
- `target::String`
- `score::Float32`
- `offset::Int`
- `orientation::String`
- `metric::String`
- `n_sites::Int`
- `p_value::Union{Nothing,Float64}`
- `adj_p_value::Union{Nothing,Float64}`
- `e_value::Union{Nothing,Float64}`
- `null_id::Union{Nothing,String}`
- `null_n::Union{Nothing,Int}`
- `null_estimator::Union{Nothing,String}`
"""
struct AnnotatedResult
    query::String
    target::String
    score::Float32
    offset::Int
    orientation::String
    metric::String
    n_sites::Int
    p_value::Union{Nothing,Float64}
    adj_p_value::Union{Nothing,Float64}
    e_value::Union{Nothing,Float64}
    null_id::Union{Nothing,String}
    null_n::Union{Nothing,Int}
    null_estimator::Union{Nothing,String}
end

function AnnotatedResult(
    result::ComparisonResult;
    p_value=nothing,
    adj_p_value=nothing,
    e_value=nothing,
    null_id=nothing,
    null_n=nothing,
    null_estimator=nothing,
)
    return AnnotatedResult(
        result.query,
        result.target,
        result.score,
        result.offset,
        result.orientation,
        result.metric,
        result.n_sites,
        p_value,
        adj_p_value,
        e_value,
        null_id,
        null_n,
        null_estimator,
    )
end

"""
    annotate_results(results, dist; effective_number_of_targets=nothing)

Annotate comparison results with p-value, adjusted p-value, and E-value from
a fitted null distribution.

Returns a vector of [`AnnotatedResult`](@ref).
"""
function annotate_results(
    results::AbstractVector{ComparisonResult},
    dist::NullDistribution;
    effective_number_of_targets::Union{Nothing,Int}=nothing,
)
    fit_result = dist.fit
    if fit_result isa GEVFitFailure
        throw(ArgumentError("Cannot annotate with a failed GEV fit: $(fit_result.message)"))
    end

    gev = fit_result::GEVFit
    n_null = dist.n_null
    effective = if effective_number_of_targets === nothing
        length(results)
    else
        effective_number_of_targets
    end
    effective >= 0 ||
        throw(ArgumentError("effective_number_of_targets must be non-negative."))

    pvalues = Vector{Float64}(undef, length(results))

    for (idx, result) in enumerate(results)
        pvalues[idx] = survival(gev, result.score)
    end

    adj = adjusted_pvalues(pvalues)

    null_id = _null_id(dist)

    annotated = Vector{AnnotatedResult}(undef, length(results))
    for idx in eachindex(results)
        r = results[idx]
        annotated[idx] = AnnotatedResult(
            r;
            p_value=pvalues[idx],
            adj_p_value=adj[idx],
            e_value=evalue(pvalues[idx], effective),
            null_id=null_id,
            null_n=n_null,
            null_estimator="genextreme",
        )
    end

    return annotated
end

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

function _null_id(dist::NullDistribution)
    parts = [
        "format_version=$(NULL_FORMAT_VERSION)",
        "strategy=$(dist.strategy)",
        "metric=$(dist.metric)",
        "n_null=$(dist.n_null)",
        "raw=$(dist.contract.raw_scores_fingerprint)",
        "seq=$(dist.contract.sequence_fingerprint)",
        "bg=$(dist.contract.background_fingerprint)",
        "contract=$(join(string.((dist.contract.search_range,
            dist.contract.window_radius, dist.contract.realign_window,
            dist.contract.min_logfpr, dist.contract.normalization_version,
            dist.contract.alignment_version)), ":"))",
    ]
    if dist.model_collection_fingerprint !== nothing
        push!(parts, "mcf=$(dist.model_collection_fingerprint)")
    end
    if dist.relation_fingerprint !== nothing
        push!(parts, "rf=$(dist.relation_fingerprint)")
    end
    dist.fit isa GEVFit && append!(
        parts,
        [
            "shape=$(dist.fit.shape)",
            "location=$(dist.fit.location)",
            "scale=$(dist.fit.scale)",
            "iterations=$(dist.fit.iterations)",
            "loglikelihood=$(dist.fit.loglikelihood)",
        ],
    )
    return bytes2hex(SHA.sha256(Vector{UInt8}(codeunits(join(parts, "|")))))
end
