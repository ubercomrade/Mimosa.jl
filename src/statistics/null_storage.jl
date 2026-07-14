# Portable null-distribution storage: TOML manifest + NPY binary blobs.
#
# Implements the profile-only v2 schema from ADR 0003. The manifest is a TOML file
# (stdlib, no external dependency) with metadata and array references.
# Binary data uses the standard NPY format (little-endian Float64).
#
# Directory bundle structure:
#   path/
#   ├── manifest.toml
#   └── data/
#       └── raw_null_scores.npy

"""
    NULL_FORMAT_VERSION

Current version of the portable null-distribution bundle format.
Value is `3`. Bundles with a different version are rejected on load.
"""
const NULL_FORMAT_VERSION = 3

"""
    savenull(path, dist::NullDistribution)

Save a null distribution to a portable directory bundle at `path`.

Creates a `manifest.toml` and a `data/raw_null_scores.npy` file in a complete
sibling staging directory, then commits the directory with an atomic rename.
"""
function savenull(path::AbstractString, dist::NullDistribution)
    dist.n_null == length(dist.raw_scores) ||
        throw(InvariantError("null distribution n_null does not match raw_scores length."))
    length(dist.pairs) == dist.n_null ||
        throw(InvariantError("null distribution pairs do not match n_null."))
    dist.n_queries >= 0 ||
        throw(InvariantError("null distribution n_queries must be non-negative."))
    isempty(dist.strategy) &&
        throw(InvariantError("null distribution strategy must not be empty."))
    isempty(dist.metric) &&
        throw(InvariantError("null distribution metric must not be empty."))
    dist.strategy == "profile" ||
        throw(InvariantError("only profile null distributions are supported."))
    dist.metric in ("co", "co_rowwise", "dice", "dice_rowwise", "cosine") ||
        throw(InvariantError("unsupported profile metric '$(dist.metric)'."))
    all(isfinite, dist.raw_scores) ||
        throw(InvariantError("null distribution raw_scores contain non-finite values."))
    _bundle_shape_payload_bytes(
        [length(dist.raw_scores)], "<f8", String(path), "raw_null_scores array"
    )

    fit_result = dist.fit
    if fit_result isa GEVFit
        all(isfinite, (fit_result.shape, fit_result.location, fit_result.scale)) ||
            throw(InvariantError("GEV fit parameters must be finite."))
        fit_result.scale > 0 || throw(InvariantError("GEV fit scale must be positive."))
        fit_result.iterations >= 0 ||
            throw(InvariantError("GEV fit iterations must be non-negative."))
        gev_params = [fit_result.shape, fit_result.location, fit_result.scale]
        estimator_type = "genextreme"
        converged = fit_result.converged
        iterations = fit_result.iterations
        loglikelihood = fit_result.loglikelihood
    else
        fit_result.iterations >= 0 ||
            throw(InvariantError("GEV failure iterations must be non-negative."))
        gev_params = [0.0, 0.0, 0.0]
        estimator_type = "failed"
        converged = false
        iterations = fit_result.iterations
        loglikelihood = 0.0
    end

    return _with_bundle_write(path) do target, stage
        npy_path = joinpath(stage, BUNDLE_DATA_DIR, "raw_null_scores.npy")
        _write_npy(npy_path, dist.raw_scores)
        checksum = _file_sha256(npy_path)
        manifest = Dict{String,Any}(
            "format" => "mimosa",
            "format_version" => NULL_FORMAT_VERSION,
            "kind" => "null_distribution",
            "strategy" => dist.strategy,
            "metric" => dist.metric,
            "estimator_type" => estimator_type,
            "genextreme_params" => gev_params,
            "genextreme_converged" => converged,
            "genextreme_iterations" => iterations,
            "genextreme_loglikelihood" => loglikelihood,
            "n_null" => dist.n_null,
            "n_queries" => dist.n_queries,
            "pairs" => [
                Dict("query" => pair.query, "target" => pair.target, "score" => pair.score) for pair in dist.pairs
            ],
            "skipped" =>
                [Dict("query" => s.query, "reason" => s.reason) for s in dist.skipped],
            "compatibility" => Dict{String,Any}(
                "format_version" => NULL_FORMAT_VERSION,
                "strategy" => dist.strategy,
                "metric" => dist.metric,
                "sequence_fingerprint" => dist.sequence_fingerprint,
                "background_fingerprint" => dist.background_fingerprint,
                "model_collection_fingerprint" =>
                    if dist.model_collection_fingerprint === nothing
                        "none"
                    else
                        dist.model_collection_fingerprint
                    end,
                "relation_fingerprint" => if dist.relation_fingerprint === nothing
                    "none"
                else
                    dist.relation_fingerprint
                end,
                "search_range" => dist.contract.search_range,
                "window_radius" => dist.contract.window_radius,
                "realign_window" => dist.contract.realign_window,
                "min_logfpr" => dist.contract.min_logfpr,
                "normalization_version" => dist.contract.normalization_version,
                "alignment_version" => dist.contract.alignment_version,
                "raw_scores_fingerprint" => dist.contract.raw_scores_fingerprint,
            ),
            "arrays" => Dict{String,Any}(
                "raw_null_scores" => Dict{String,Any}(
                    "file" => "data/raw_null_scores.npy",
                    "dtype" => "<f8",
                    "shape" => [length(dist.raw_scores)],
                    "checksum" => "sha256:$checksum",
                ),
            ),
        )
        _write_bundle_manifest(joinpath(stage, BUNDLE_MANIFEST_NAME), manifest)
        return target
    end
end

"""
    loadnull(path)

Load a null distribution from a directory bundle at `path`.

Validates the manifest format version and checksums. Returns a
[`NullDistribution`](@ref).
"""
function loadnull(path::AbstractString)
    try
        manifest = _read_bundle_manifest(
            path, NULL_FORMAT_VERSION; expected_kind="null_distribution"
        )
        strategy = _required_manifest_string(manifest, "strategy", path, "null manifest")
        strategy == "profile" ||
            throw(_bundle_error(path, "unsupported null strategy '$strategy'."))
        metric = _required_manifest_string(manifest, "metric", path, "null manifest")
        n_null = _required_manifest_int(
            manifest,
            "n_null",
            path,
            "null manifest";
            minimum=0,
            maximum=MAX_BUNDLE_ELEMENTS,
        )
        n_queries = _required_manifest_int(
            manifest,
            "n_queries",
            path,
            "null manifest";
            minimum=0,
            maximum=MAX_BUNDLE_ELEMENTS,
        )

        arrays = _required_manifest_table(manifest, "arrays", path, "null manifest")
        spec = _parse_bundle_array(path, arrays, "raw_null_scores", path)
        spec.dtype == "<f8" ||
            throw(_bundle_error(path, "raw_null_scores must use dtype '<f8'."))
        spec.shape == [n_null] ||
            throw(_bundle_error(path, "raw_null_scores shape does not match n_null."))
        npy_path = _validate_bundle_array_checksum(path, spec, path)
        raw_scores = _read_npy_f64(npy_path; expected_shape=spec.shape)
        all(isfinite, raw_scores) ||
            throw(_bundle_error(path, "raw_null_scores contains non-finite values."))

        gev_params = _required_manifest_floats(
            manifest, "genextreme_params", path, "null manifest"; expected_length=3
        )
        estimator_type = _required_manifest_string(
            manifest, "estimator_type", path, "null manifest"
        )
        converged = _required_manifest_bool(
            manifest, "genextreme_converged", path, "null manifest"
        )
        iterations = _required_manifest_int(
            manifest,
            "genextreme_iterations",
            path,
            "null manifest";
            minimum=0,
            maximum=MAX_BUNDLE_ELEMENTS,
        )
        loglikelihood = _required_manifest_float(
            manifest, "genextreme_loglikelihood", path, "null manifest"
        )

        fit_result = if estimator_type == "genextreme"
            gev_params[3] > 0 || throw(_bundle_error(path, "GEV scale must be positive."))
            GEVFit(
                gev_params[1],
                gev_params[2],
                gev_params[3],
                converged,
                iterations,
                loglikelihood,
            )
        elseif estimator_type == "failed"
            GEVFitFailure(
                "Stored null distribution has a failed GEV fit.",
                length(raw_scores),
                iterations,
            )
        else
            throw(_bundle_error(path, "unsupported estimator_type '$estimator_type'."))
        end

        skipped_raw = get(manifest, "skipped", nothing)
        skipped_raw isa AbstractVector ||
            throw(_bundle_error(path, "null manifest 'skipped' must be an array."))
        skipped = NamedTuple{(:query, :reason),Tuple{String,String}}[]
        for (index, item) in enumerate(skipped_raw)
            item isa AbstractDict ||
                throw(_bundle_error(path, "skipped entry $index must be a TOML table."))
            query = _required_manifest_string(item, "query", path, "skipped entry $index")
            reason = _required_manifest_string(item, "reason", path, "skipped entry $index")
            push!(skipped, (query=query, reason=reason))
        end

        pairs_raw = get(manifest, "pairs", nothing)
        pairs_raw isa AbstractVector ||
            throw(_bundle_error(path, "null manifest 'pairs' must be an array."))
        pairs = NullPair[]
        for (index, item) in enumerate(pairs_raw)
            item isa AbstractDict ||
                throw(_bundle_error(path, "pair entry $index must be a TOML table."))
            query = _required_manifest_string(item, "query", path, "pair entry $index")
            target = _required_manifest_string(item, "target", path, "pair entry $index")
            score = _required_manifest_float(item, "score", path, "pair entry $index")
            push!(pairs, NullPair(query, target, score))
        end
        length(pairs) == n_null ||
            throw(_bundle_error(path, "pair count does not match n_null."))

        compat = _required_manifest_table(manifest, "compatibility", path, "null manifest")
        compat_version = _required_manifest_int(
            compat,
            "format_version",
            path,
            "compatibility metadata";
            minimum=1,
            maximum=NULL_FORMAT_VERSION,
        )
        compat_version == NULL_FORMAT_VERSION ||
            throw(_bundle_error(path, "unsupported compatibility metadata version."))
        compat_strategy = _required_manifest_string(
            compat, "strategy", path, "compatibility metadata"
        )
        compat_metric = _required_manifest_string(
            compat, "metric", path, "compatibility metadata"
        )
        compat_strategy == strategy || throw(
            _bundle_error(path, "compatibility strategy disagrees with null manifest.")
        )
        compat_metric == metric ||
            throw(_bundle_error(path, "compatibility metric disagrees with null manifest."))
        seq_fp = _required_manifest_string(
            compat, "sequence_fingerprint", path, "compatibility metadata"
        )
        bg_fp = _required_manifest_string(
            compat, "background_fingerprint", path, "compatibility metadata"
        )
        mcf = _required_manifest_string(
            compat, "model_collection_fingerprint", path, "compatibility metadata"
        )
        rf = _required_manifest_string(
            compat, "relation_fingerprint", path, "compatibility metadata"
        )
        search_range = _required_manifest_int(
            compat, "search_range", path, "compatibility metadata"; minimum=0
        )
        window_radius = _required_manifest_int(
            compat, "window_radius", path, "compatibility metadata"; minimum=0
        )
        realign_window = _required_manifest_int(
            compat, "realign_window", path, "compatibility metadata"; minimum=0
        )
        min_logfpr = _required_manifest_float(
            compat, "min_logfpr", path, "compatibility metadata"
        )
        normalization_version = _required_manifest_string(
            compat, "normalization_version", path, "compatibility metadata"
        )
        alignment_version = _required_manifest_string(
            compat, "alignment_version", path, "compatibility metadata"
        )
        raw_scores_fingerprint = _required_manifest_string(
            compat, "raw_scores_fingerprint", path, "compatibility metadata"
        )
        actual_raw_scores_fingerprint = content_fingerprint(raw_scores)
        raw_scores_fingerprint == actual_raw_scores_fingerprint ||
            throw(_bundle_error(path, "raw score fingerprint does not match payload."))
        contract = ProfileComparisonContract(
            metric,
            search_range,
            window_radius,
            realign_window,
            Float32(min_logfpr),
            normalization_version,
            alignment_version,
            seq_fp,
            bg_fp,
            raw_scores_fingerprint,
        )

        return NullDistribution(
            strategy,
            metric,
            fit_result,
            raw_scores,
            pairs,
            n_null,
            n_queries,
            skipped,
            mcf != "none" ? mcf : nothing,
            rf != "none" ? rf : nothing,
            seq_fp,
            bg_fp,
            contract,
        )
    catch err
        _rethrow_bundle_error(path, err)
    end
end
