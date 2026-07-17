"""
    Mimosa

A Julia package for motif scanning, comparison, and statistical evaluation.

Supports six model families (PWM, PFM, BaMM, SiteGA, Dimont, Slim) with
profile-based comparison, site extraction, PFM reconstruction,
native GEV null distributions, BH FDR, E-values, portable storage bundles,
content-based cache, serial and threaded parallelism, and a thin CLI adapter.

See `docs/src/api.md` and `docs/src/data_layout.md` for usage and coordinate conventions.
"""
module Mimosa

include("errors.jl")
include("parallel/parallel.jl")
include("models/models.jl")
include("sequences/sequences.jl")
include("scanning/scanning.jl")
include("profiles/profiles.jl")
include("io/io.jl")
include("comparison/comparison.jl")
include("sites/sites.jl")
include("statistics/statistics.jl")
include("cache/cache.jl")
include("serialization.jl")
include("cli.jl")
include("precompile.jl")

"""
    main(args=ARGS)

Run the Mimosa command-line interface and return its exit code.
"""
function (@main)(args::Vector{String}=ARGS)
    return cli_main(args)
end

export readmodel,
    read_scores,
    read_bamm,
    read_sitega,
    read_dimont,
    read_slim,
    write_sitega,
    writemodel,
    readsequences,
    compare,
    to_json,
    to_dict
export MimosaError,
    ModelFormatError, ModelDimensionError, InvariantError, ModelInterfaceError
export AbstractProfileSource,
    AbstractMotifModel,
    PWM,
    BaMM,
    SiteGA,
    Dimont,
    Slim,
    pcm_to_pfm,
    pfm_to_pwm,
    pwm_from_pfm,
    extend_pwm_with_n,
    modelname,
    left_context,
    right_context,
    scan_kernel!,
    site_start_offset,
    validate_model,
    ComparisonResult,
    metric_name,
    read_fasta

# Sequence / scanning exports
export EncodedSequenceBatch,
    nsequences,
    seqlength,
    sequence,
    empty_sequence_batch,
    make_random_sequences,
    encode_base,
    encode_sequence,
    reverse_complement,
    reverse_complement!,
    to_padded,
    from_padded,
    N_CODE
export RaggedArray, nrows, rowlength, row, build_ragged, empty_ragged
export StrandPolicy, ForwardOnly, ReverseOnly, BestStrand, BothStrands, StrandPair
export scan,
    scan!,
    scan_forward!,
    scan_reverse!,
    scan_best_strand!,
    best_hits!,
    scan_both!,
    npositions,
    scan_result_lengths

# Model exports
export ScoreProfile,
    scorebounds, profile_bundle, motif_length, window_size, scorematrix, scoretype
export is_scannable

# Sites and PFM reconstruction exports
export SiteSelector,
    BestPerSequence,
    ThresholdHits,
    TopFractionHits,
    SiteHit,
    SiteCollection,
    selectsites,
    reconstruct_pfm,
    extract_site_matrix,
    build_pcm,
    site_strings,
    sort_hits!,
    select_top_fraction,
    empty_site_collection

# Profile comparison exports
export AbstractProfileMetric,
    OverlapCoefficient,
    OverlapCoefficientRowwise,
    DiceSimilarity,
    DiceSimilarityRowwise,
    CosineSimilarityProfile,
    parse_profile_metric,
    ProfileConfig,
    profile_compare,
    PreparedProfile,
    prepare_profile,
    LogTailTable,
    EmpiricalLogTail,
    fit,
    lookup_score,
    transform_scores,
    flatten_bundle,
    normalize_bundle,
    AnchorCSR,
    build_anchor_csr,
    collect_best_anchors,
    collect_threshold_anchors,
    collect_anchors,
    score_shift

# Statistics exports
export GEVFit,
    GEVFitFailure,
    GEVFitResult,
    fit_gev,
    survival,
    cdf,
    scipy_params,
    BenjaminiHochberg,
    adjusted_pvalues,
    evalue,
    pvalue,
    GroupRelations,
    parse_group_relations,
    eligible_targets,
    NullDistribution,
    NullPair,
    ProfileComparisonContract,
    NullBuildConfig,
    NullBuildResult,
    AnnotatedResult,
    build_null,
    annotate_results,
    ANNOTATED_RESULT_SCHEMA_VERSION,
    savenull,
    loadnull,
    NULL_FORMAT_VERSION

# Parallelism exports
export ExecutionPolicy, SerialExecution, ThreadedExecution

# Cache exports
export Cache,
    cache_key,
    prepared_profile_cache_key,
    cache_has,
    cache_get,
    cache_get_meta,
    cache_set,
    clearcache,
    content_fingerprint,
    model_fingerprint,
    model_collection_fingerprint,
    sequence_fingerprint,
    MODEL_FORMAT_VERSION

end # module
