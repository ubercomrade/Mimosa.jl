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
- `n_models::Int`: number of source models.
- `model_type::String`: model family used to build the null.
- `shuffle::Bool`: whether PWM shuffling was enabled.
- `seed::Int`: seed used for pair sampling and shuffling.
- `sampling_version::String`: random sampling algorithm version.
- `model_collection_fingerprint::Union{Nothing,String}`: SHA-256 of model collection.
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
    n_models::Int
    model_type::String
    shuffle::Bool
    seed::Int
    sampling_version::String
    model_collection_fingerprint::Union{Nothing,String}
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
    n_models::Int,
    model_type::String,
    shuffle::Bool,
    seed::Int,
    sampling_version::String,
    model_collection_fingerprint,
    sequence_fingerprint::String,
    background_fingerprint::String,
)
    n_null >= 0 || throw(InvariantError("null distribution n_null must be non-negative."))
    n_null == length(raw_scores) ||
        throw(InvariantError("null distribution n_null does not match raw_scores length."))
    n_models >= 2 ||
        throw(InvariantError("null distribution requires at least two source models."))
    length(pairs) == n_null ||
        throw(InvariantError("null distribution pairs do not match n_null."))
    isempty(model_type) &&
        throw(InvariantError("null distribution model_type must not be empty."))
    seed >= 0 || throw(InvariantError("null distribution seed must be non-negative."))
    isempty(sampling_version) &&
        throw(InvariantError("null distribution sampling_version must not be empty."))
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
        n_models,
        model_type,
        shuffle,
        seed,
        sampling_version,
        model_collection_fingerprint,
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
- `n_samples::Int`: number of random pair comparisons.
- `shuffle::Bool`: whether PWM models are shuffled before each comparison.
- `seed::Int`: random seed for pair selection and PWM shuffling.
"""
struct NullBuildConfig{M<:AbstractProfileMetric}
    metric::M
    n_samples::Int
    shuffle::Bool
    seed::Int

    function NullBuildConfig(
        metric::M, n_samples::Int, shuffle::Bool, seed::Int
    ) where {M<:AbstractProfileMetric}
        n_samples > 0 || throw(ArgumentError("n_samples must be positive."))
        seed >= 0 || throw(ArgumentError("seed must be non-negative."))
        return new{M}(metric, n_samples, shuffle, seed)
    end
end

function NullBuildConfig(;
    metric=nothing, n_samples::Int=2000, shuffle::Bool=false, seed::Int=127
)
    resolved_metric = _resolve_profile_metric(isnothing(metric) ? :co : metric)
    return NullBuildConfig(resolved_metric, n_samples, shuffle, seed)
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
    build_null(models; sequences=batch, metric=:co, n_samples=2000,
               shuffle=false, seed=127, outer_execution=SerialExecution(),
               scan_execution=SerialExecution())

Build a null distribution from randomly sampled query-target comparisons.

Each sample selects two distinct models uniformly, optionally shuffles both
models when they are PWMs, and collects their comparison score. Sampling is
with replacement between iterations. The pooled scores are fitted to a GEV
distribution.

# Arguments
- `models::AbstractVector`: vector of motif models (e.g. `PWM`).
- `metric`: comparison metric.
- `n_samples`: number of random comparisons (default 2000).
- `shuffle`: shuffle PWM columns and A/C/G/T weights within each column.
- `seed`: random seed for reproducible sampling and shuffling.
- `outer_execution`: [`ExecutionPolicy`](@ref) for parallel comparison of query-target
  pairs. Default `SerialExecution()`.
- `scan_execution`: policy for scanning sequences and applying normalization
  within one pair. It must not be multi-threaded together with `outer_execution`.

Returns a [`NullBuildResult`](@ref).

Under `ThreadedExecution`, comparisons are processed in parallel at the
top level. Results are collected into pre-allocated slots indexed by
the original comparison order, so the pooled score order and fit are
identical to `SerialExecution`.
"""
function build_null(
    models::AbstractVector;
    metric=nothing,
    n_samples::Int=2000,
    shuffle::Bool=false,
    seed::Int=127,
    outer_execution::ExecutionPolicy=SerialExecution(),
    scan_execution::ExecutionPolicy=SerialExecution(),
    sequences::EncodedSequenceBatch,
    background::Union{Nothing,EncodedSequenceBatch}=nothing,
    search_range::Int=10,
    window_radius::Int=10,
    realign_window::Int=3,
    min_logfpr::Real=0.0,
    normalization::AbstractNormalizationStrategy=HybridEmpiricalLogTail(),
)
    config = NullBuildConfig(;
        metric=metric, n_samples=n_samples, shuffle=shuffle, seed=seed
    )
    return build_null(
        models,
        config;
        outer_execution=outer_execution,
        scan_execution=scan_execution,
        sequences=sequences,
        background=background,
        search_range=search_range,
        window_radius=window_radius,
        realign_window=realign_window,
        min_logfpr=min_logfpr,
        normalization=normalization,
    )
end

function build_null(
    models::AbstractVector,
    config::NullBuildConfig{<:AbstractProfileMetric};
    outer_execution::ExecutionPolicy=SerialExecution(),
    scan_execution::ExecutionPolicy=SerialExecution(),
    sequences::EncodedSequenceBatch,
    background::Union{Nothing,EncodedSequenceBatch}=nothing,
    search_range::Int=10,
    window_radius::Int=10,
    realign_window::Int=3,
    min_logfpr::Real=0.0,
    normalization::AbstractNormalizationStrategy=HybridEmpiricalLogTail(),
)
    _validate_execution_levels(outer_execution, scan_execution)
    length(models) >= 2 ||
        throw(ArgumentError("at least two models are required for null construction."))
    all(model -> model isa AbstractMotifModel, models) ||
        throw(ArgumentError("models must contain only AbstractMotifModel values."))
    for model in models
        validate_model(model; capability=:cache)
    end
    names = String[modelname(model) for model in models]
    length(unique(names)) == length(names) ||
        throw(ArgumentError("model names must be unique for null construction."))
    all(!isempty, names) || throw(ArgumentError("model names must not be empty."))
    model_types = unique(_null_model_type(model) for model in models)
    length(model_types) == 1 ||
        throw(ArgumentError("all models must belong to the same model family."))
    search_range >= 0 || throw(ArgumentError("search_range must be non-negative."))
    window_radius >= 0 || throw(ArgumentError("window_radius must be non-negative."))
    realign_window >= 0 || throw(ArgumentError("realign_window must be non-negative."))
    isfinite(min_logfpr) || throw(ArgumentError("min_logfpr must be finite."))
    prepared = Vector{Union{Nothing,PreparedProfile}}(undef, length(models))
    fill!(prepared, nothing)
    reusable = findall(model -> !config.shuffle || !(model isa PWM), models)
    _parallel_for(outer_execution, length(reusable)) do i
        model_index = reusable[i]
        prepared[model_index] = prepare_profile(
            models[model_index],
            sequences;
            background=background,
            min_logfpr=min_logfpr,
            normalization=normalization,
            scan_execution=scan_execution,
        )
        return nothing
    end

    work_items = _null_work_items(length(models), config)
    raw_scores = Vector{Float64}(undef, config.n_samples)
    pairs = Vector{NullPair}(undef, config.n_samples)
    _parallel_for(outer_execution, config.n_samples) do i
        item = work_items[i]
        query = models[item.query]
        target = models[item.target]
        query_profile = if isnothing(prepared[item.query])
            prepare_profile(
                _shuffle_null_model(query, item.query_seed),
                sequences;
                background=background,
                min_logfpr=min_logfpr,
                normalization=normalization,
                scan_execution=scan_execution,
            )
        else
            prepared[item.query]::PreparedProfile
        end
        target_profile = if isnothing(prepared[item.target])
            prepare_profile(
                _shuffle_null_model(target, item.target_seed),
                sequences;
                background=background,
                min_logfpr=min_logfpr,
                normalization=normalization,
                scan_execution=scan_execution,
            )
        else
            prepared[item.target]::PreparedProfile
        end
        result = compare(
            query_profile,
            target_profile;
            metric=config.metric,
            search_range=search_range,
            window_radius=window_radius,
            realign_window=realign_window,
        )
        score = Float64(result.score)
        raw_scores[i] = score
        return pairs[i] = NullPair(
            String(modelname(query)), String(modelname(target)), score
        )
    end

    fit_result = fit_gev(raw_scores)
    seq_fingerprint = sequence_fingerprint(sequences)
    bg_fingerprint = isnothing(background) ? "none" : sequence_fingerprint(background)
    metric = metric_name(config.metric)
    dist = NullDistribution(
        "profile",
        metric,
        fit_result,
        raw_scores,
        pairs,
        length(raw_scores),
        length(models),
        only(model_types),
        config.shuffle,
        config.seed,
        "random-ordered-pairs-v1",
        model_collection_fingerprint(AbstractProfileSource[models...]),
        seq_fingerprint,
        bg_fingerprint,
        ProfileComparisonContract(
            metric,
            search_range,
            window_radius,
            realign_window,
            Float32(min_logfpr),
            normalization_fingerprint(normalization),
            "profile-alignment-v1",
            seq_fingerprint,
            bg_fingerprint,
            content_fingerprint(raw_scores),
        ),
    )
    return NullBuildResult(dist, config.n_samples)
end

struct _NullWorkItem
    query::Int
    target::Int
    query_seed::UInt64
    target_seed::UInt64
end

function _null_work_items(n_models::Int, config::NullBuildConfig)
    rng = Random.MersenneTwister(config.seed)
    items = Vector{_NullWorkItem}(undef, config.n_samples)
    for i in eachindex(items)
        query = rand(rng, 1:n_models)
        target = rand(rng, 1:(n_models - 1))
        target >= query && (target += 1)
        items[i] = _NullWorkItem(query, target, rand(rng, UInt64), rand(rng, UInt64))
    end
    return items
end

_null_model_type(model::AbstractMotifModel) = lowercase(String(nameof(typeof(model))))

_shuffle_null_model(model::AbstractMotifModel, ::UInt64) = model

function _shuffle_null_model(model::PWM, seed::UInt64)
    rng = Random.MersenneTwister(seed)
    column_order = randperm(rng, motif_length(model))
    representation = copy(model.representation[:, column_order])
    for column in axes(representation, 2)
        base_order = randperm(rng, NUCLEOTIDE_CARDINALITY)
        representation[1:NUCLEOTIDE_CARDINALITY, column] = representation[
            base_order, column
        ]
        representation[5, column] = minimum(
            @view representation[1:NUCLEOTIDE_CARDINALITY, column]
        )
    end
    return PWM(model.name, representation, model.background)
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
        result.metric == dist.metric || throw(
            ArgumentError(
                "result metric '$(result.metric)' does not match null metric '$(dist.metric)'.",
            ),
        )
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
        "n_models=$(dist.n_models)",
        "model_type=$(dist.model_type)",
        "shuffle=$(dist.shuffle)",
        "seed=$(dist.seed)",
        "sampling=$(dist.sampling_version)",
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
