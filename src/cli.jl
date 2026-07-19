# CLI: thin adapter from command-line arguments to the public Mimosa API.
#
# ArgParse owns command-line syntax and help generation. Each subcommand maps
# to a typed configuration, calls the public API, and
# serializes results as JSON to stdout. Diagnostics and progress go to stderr.
#
# Exit codes:
#   0 = success
#   1 = usage error / invalid arguments
#   2 = runtime error (file not found, malformed input, etc.)
#
# JSON output goes to stdout only. Logs and errors go to stderr only.

using ArgParse

const CLI_VERSION = string(Base.pkgversion(@__MODULE__))

const MODEL_TYPES = ["pwm", "bamm", "sitega", "dimont", "slim"]
const PROFILE_MODEL_TYPES = ["scores", MODEL_TYPES...]
const PROFILE_METRICS = ["co", "co_rowwise", "dice", "dice_rowwise", "cosine"]

struct CLIError <: Exception
    message::String
end

struct CLIParsed
    command::String
    positional::Vector{String}
    options::Dict{String,String}
    flags::Set{String}
end

function _cli_int(parsed::CLIParsed, name::String, default::String; minimum=nothing)
    value = tryparse(Int, get(parsed.options, name, default))
    value === nothing && throw(CLIError("--$name must be an integer."))
    minimum !== nothing &&
        value < minimum &&
        throw(CLIError("--$name must be at least $minimum."))
    return value
end

function _cli_float32(parsed::CLIParsed, name::String, default::String)
    value = tryparse(Float32, get(parsed.options, name, default))
    value === nothing && throw(CLIError("--$name must be a finite Float32."))
    isfinite(value) || throw(CLIError("--$name must be a finite Float32."))
    return value
end

function CLIParsed(command::String)
    return CLIParsed(command, String[], Dict{String,String}(), Set{String}())
end

# ── Type and model resolution ────────────────────────────────────────────────
#
# Shared helpers used by command runners to read models and resolve sequences.
# These bridge the parsed CLI strings to typed Mimosa API objects.

"""
    _read_typed_model(path, model_type; kwargs...)

Read a model file using the specified type string.
"""
function _read_typed_model(path::AbstractString, model_type::AbstractString; kwargs...)
    if model_type == "scores"
        return read_scores(path)
    elseif model_type == "pwm"
        return readmodel(path; kwargs...)
    elseif model_type == "bamm"
        order_val = get(kwargs, :order, nothing)
        return read_bamm(path; order=order_val)
    elseif model_type == "sitega"
        return read_sitega(path)
    elseif model_type == "dimont"
        return read_dimont(path)
    elseif model_type == "slim"
        return read_slim(path)
    else
        throw(CLIError("unknown model type: $(model_type)"))
    end
end

"""
    _resolve_sequences(fasta_path, num_sequences, seq_length, seed)

Resolve sequences from a FASTA file or generate random ones.
"""
function _resolve_sequences(
    fasta_path::Union{AbstractString,Nothing},
    num_sequences::Int,
    seq_length::Int,
    seed::Int,
)
    if fasta_path !== nothing
        batch, _ = read_fasta(fasta_path)
        return batch
    end
    return make_random_sequences(num_sequences, seq_length; seed=seed)
end

function _execution_policy(parsed::CLIParsed)
    threads_str = get(parsed.options, "threads", get(parsed.options, "jobs", "1"))
    requested = tryparse(Int, threads_str)
    requested === nothing && throw(CLIError("--threads must be a positive integer."))
    requested > 0 || throw(CLIError("--threads must be a positive integer."))

    available = Threads.nthreads()
    requested <= available || throw(
        CLIError(
            "--threads=$requested exceeds the $available Julia thread(s) available; " *
            "start Julia with --threads=$requested or set JULIA_NUM_THREADS=$requested.",
        ),
    )
    return requested == 1 ? SerialExecution() : ThreadedExecution(requested)
end

# ── JSON output ─────────────────────────────────────────────────────────────

function _json_value(v::Nothing)
    return "null"
end

function _json_value(v::AbstractString)
    return _json_string(v)
end

function _json_value(v::AbstractFloat)
    return _json_float(Float64(v))
end

function _json_value(v::Real)
    return string(v)
end

function _json_value(v::Bool)
    return v ? "true" : "false"
end

function _json_value(v::AbstractVector)
    if isempty(v)
        return "[]"
    end
    parts = [_json_value(x) for x in v]
    return "[" * join(parts, ", ") * "]"
end

function _json_dict(d::Dict{String})
    keys_sorted = sort!(collect(keys(d)))
    parts = String[]
    for k in keys_sorted
        push!(parts, _json_string(k) * ": " * _json_value(d[k]))
    end
    return "{" * join(parts, ", ") * "}"
end

function _println_json(d::Dict{String})
    println(stdout, _json_dict(d))
    return nothing
end

function _validate_null_compatibility(
    dist::NullDistribution;
    strategy::AbstractString,
    metric::AbstractString,
    sequences::Union{Nothing,EncodedSequenceBatch}=nothing,
    background::Union{Nothing,EncodedSequenceBatch}=nothing,
    search_range::Int=10,
    window_radius::Int=10,
    realign_window::Int=3,
    min_logfpr::Float32=0.0f0,
    model_types::Union{Nothing,Tuple{String,String}}=nothing,
)
    dist.strategy == strategy || throw(
        CLIError(
            "null distribution strategy '$(dist.strategy)' is incompatible with $strategy comparison.",
        ),
    )
    dist.metric == metric || throw(
        CLIError(
            "null distribution metric '$(dist.metric)' is incompatible with requested metric '$metric'.",
        ),
    )
    contract = dist.contract
    contract.search_range == search_range || throw(
        CLIError("null distribution search range is incompatible with this comparison.")
    )
    contract.window_radius == window_radius || throw(
        CLIError("null distribution window radius is incompatible with this comparison."),
    )
    contract.realign_window == realign_window || throw(
        CLIError(
            "null distribution realignment window is incompatible with this comparison."
        ),
    )
    contract.min_logfpr == min_logfpr || throw(
        CLIError("null distribution minimum log-FPR is incompatible with this comparison."),
    )

    expected_sequences = isnothing(sequences) ? "none" : sequence_fingerprint(sequences)
    expected_background = isnothing(background) ? "none" : sequence_fingerprint(background)
    dist.sequence_fingerprint == expected_sequences || throw(
        CLIError(
            "null distribution sequence fingerprint is incompatible with this comparison.",
        ),
    )
    dist.background_fingerprint == expected_background || throw(
        CLIError(
            "null distribution background fingerprint is incompatible with this comparison.",
        ),
    )
    if model_types !== nothing
        all(==(dist.model_type), model_types) || throw(
            CLIError(
                "null distribution model type '$(dist.model_type)' is incompatible with compared model types '$(join(model_types, ", "))'.",
            ),
        )
    end
    return nothing
end

function _annotate_cli_result(
    result::ComparisonResult,
    parsed::CLIParsed;
    strategy::AbstractString,
    metric::AbstractString,
    sequences::Union{Nothing,EncodedSequenceBatch}=nothing,
    background::Union{Nothing,EncodedSequenceBatch}=nothing,
    search_range::Int=10,
    window_radius::Int=10,
    realign_window::Int=3,
    min_logfpr::Float32=0.0f0,
    model_types::Union{Nothing,Tuple{String,String}}=nothing,
)
    "pvalue" in parsed.flags || return result
    haskey(parsed.options, "null-distribution") ||
        throw(CLIError("--pvalue requires an explicit --null-distribution bundle."))
    dist = loadnull(parsed.options["null-distribution"])
    _validate_null_compatibility(
        dist;
        strategy=strategy,
        metric=metric,
        sequences=sequences,
        background=background,
        search_range=search_range,
        window_radius=window_radius,
        realign_window=realign_window,
        min_logfpr=min_logfpr,
        model_types=model_types,
    )
    effective = if haskey(parsed.options, "effective-number-of-targets")
        value = tryparse(Int, parsed.options["effective-number-of-targets"])
        value === nothing &&
            throw(CLIError("--effective-number-of-targets must be a positive integer."))
        value > 0 ||
            throw(CLIError("--effective-number-of-targets must be a positive integer."))
        value
    else
        nothing
    end
    return only(annotate_results([result], dist; effective_number_of_targets=effective))
end

# ── Command runners (Layer 3) ──────────────────────────────────────────────
#
# Each runner takes a parsed CLIParsed, validates scientific options,
# calls the Mimosa public API, and serializes results as JSON to stdout.
# Runners do NOT re-parse arguments; they consume the typed CLIParsed struct.

# ── Command: profile ─────────────────────────────────────────────────────────

function _print_profile_help(io::IO)
    println(
        io,
        "Usage: mimosa profile <model1> <model2> --model1-type <type> --model2-type <type> [options]",
    )
    println(io, "")
    println(
        io, "Compare motifs via score profiles: precomputed scores or profiles from scans."
    )
    println(io, "")
    println(io, "Required arguments:")
    println(io, "  model1                   Path to first model or score-profile file")
    println(io, "  model2                   Path to second model or score-profile file")
    println(io, "  --model1-type <type>     Type: $(join(PROFILE_MODEL_TYPES, ", "))")
    println(io, "  --model2-type <type>     Type: $(join(PROFILE_MODEL_TYPES, ", "))")
    println(io, "")
    println(io, "Profile comparison options:")
    println(
        io,
        "  --metric <name>          Metric: $(join(PROFILE_METRICS, ", ")) (default: co)",
    )
    println(io, "  --search-range <n>       Max site-center shift (default: 10)")
    println(
        io, "  --window-radius <n>      Window radius in profile positions (default: 10)"
    )
    println(io, "  --realign-window <n>     Local realignment half-width (default: 3)")
    println(io, "  --min-logfpr <f>         Threshold logFPR (0 = best site per sequence)")
    println(io, "")
    println(io, "Sequence options:")
    println(io, "  --fasta <path>            FASTA for motif scanning")
    println(io, "  --background <path>       FASTA for normalization calibration")
    println(io, "  --num-sequences <n>       Random sequences if no FASTA (default: 1000)")
    println(io, "  --seq-length <n>          Random sequence length (default: 200)")
    println(io, "  --seed <n>                Random seed (default: 127)")
    println(io, "  --background-freq <f>      Background freq for PWM (default: 0.25)")
    println(
        io, "  --pvalue                   Annotate using an explicit compatible null bundle"
    )
    println(
        io, "  --null-distribution <p>    Portable null-distribution bundle for --pvalue"
    )
    println(io, "  --effective-number-of-targets <n>  E-value target-count override")
    println(io, "")
    println(io, "Technical options:")
    println(
        io,
        "  --threads <n>             Worker threads to use (default: 1; runtime must provide them)",
    )
    println(
        io, "  --cache-dir <path>        Persist prepared profiles in this cache directory"
    )
    println(io, "  --quiet                   Suppress informational output")
    println(io, "  --verbose                 Verbose diagnostics to stderr")
    return nothing
end

function _run_profile(parsed::CLIParsed)
    if length(parsed.positional) != 2
        _print_profile_help(stderr)
        throw(CLIError("profile requires two positional arguments."))
    end
    path1 = parsed.positional[1]
    path2 = parsed.positional[2]

    haskey(parsed.options, "model1-type") || throw(CLIError("--model1-type is required."))
    haskey(parsed.options, "model2-type") || throw(CLIError("--model2-type is required."))

    type1 = parsed.options["model1-type"]
    type2 = parsed.options["model2-type"]
    type1 in PROFILE_MODEL_TYPES ||
        throw(CLIError("--model1-type must be one of: $(join(PROFILE_MODEL_TYPES, ", "))"))
    type2 in PROFILE_MODEL_TYPES ||
        throw(CLIError("--model2-type must be one of: $(join(PROFILE_MODEL_TYPES, ", "))"))

    metric = get(parsed.options, "metric", "co")
    metric in PROFILE_METRICS ||
        throw(CLIError("--metric must be one of: $(join(PROFILE_METRICS, ", "))"))

    search_range = _cli_int(parsed, "search-range", "10"; minimum=0)
    window_radius = _cli_int(parsed, "window-radius", "10"; minimum=0)
    realign_window = _cli_int(parsed, "realign-window", "3"; minimum=0)
    min_logfpr_str = get(parsed.options, "min-logfpr", nothing)
    min_logfpr =
        min_logfpr_str === nothing ? 0.0f0 : _cli_float32(parsed, "min-logfpr", "0.0")
    isfinite(min_logfpr) || throw(CLIError("--min-logfpr must be finite."))
    seed = _cli_int(parsed, "seed", "127")
    num_seq = _cli_int(parsed, "num-sequences", "1000"; minimum=1)
    seq_len = _cli_int(parsed, "seq-length", "200"; minimum=1)
    bg_freq = _cli_float32(parsed, "background-freq", "0.25")
    fasta = get(parsed.options, "fasta", nothing)
    bg_fasta = get(parsed.options, "background", nothing)
    execution = _execution_policy(parsed)
    cache =
        haskey(parsed.options, "cache-dir") ? Cache(parsed.options["cache-dir"]) : nothing

    model1 = _read_typed_model(path1, type1; background=bg_freq)
    model2 = _read_typed_model(path2, type2; background=bg_freq)

    # If both are ScoreProfile, compare directly.
    sequences = nothing
    bg_sequences = nothing
    if model1 isa ScoreProfile && model2 isa ScoreProfile
        result = compare(
            model1,
            model2;
            metric=metric,
            search_range=search_range,
            window_radius=window_radius,
            realign_window=realign_window,
            min_logfpr=min_logfpr,
            cache=cache,
        )
    else
        # Motif-derived profiles: scan → normalize → compare
        sequences = _resolve_sequences(fasta, num_seq, seq_len, seed)
        bg_sequences = bg_fasta !== nothing ? read_fasta(bg_fasta)[1] : nothing
        result = compare(
            model1,
            model2,
            sequences;
            metric=metric,
            search_range=search_range,
            window_radius=window_radius,
            realign_window=realign_window,
            min_logfpr=min_logfpr,
            background=bg_sequences,
            execution=execution,
            cache=cache,
        )
    end

    annotated = _annotate_cli_result(
        result,
        parsed;
        strategy="profile",
        metric=metric,
        sequences=sequences,
        background=bg_sequences,
        search_range=search_range,
        window_radius=window_radius,
        realign_window=realign_window,
        min_logfpr=min_logfpr,
        model_types=(type1, type2),
    )
    _println_json(to_dict(annotated))
    return 0
end

# ── Command: build-null ──────────────────────────────────────────────────────

function _print_build_null_help(io::IO)
    println(
        io,
        "Usage: mimosa build-null <motifs-dir> --model-type <type> --output <path> [options]",
    )
    println(io, "")
    println(io, "Build a pooled null distribution from unrelated motif comparisons.")
    println(io, "")
    println(io, "Required arguments:")
    println(io, "  motifs-dir                Directory containing motif files")
    println(io, "  --model-type <type>       Motif format: $(join(MODEL_TYPES, ", "))")
    println(io, "  --output <path>           Output path for null distribution")
    println(io, "")
    println(io, "Comparison options:")
    println(io, "  --metric <name>           Profile metric (default: co)")
    println(io, "  --fasta <path>            FASTA for profile scanning")
    println(io, "  --num-sequences <n>       Random sequences (default: 1000)")
    println(io, "  --seq-length <n>          Random sequence length (default: 200)")
    println(io, "  --seed <n>                Random seed (default: 127)")
    println(io, "  --num-samples <n>         Random comparisons (default: 2000)")
    println(io, "  --shuffle                 Shuffle PWM models before each comparison")
    println(io, "  --search-range <n>        Max shift (default: 10)")
    println(io, "  --window-radius <n>       Window radius (default: 10)")
    println(io, "  --realign-window <n>      Realignment window (default: 3)")
    println(io, "  --min-logfpr <f>          Threshold logFPR")
    println(io, "")
    println(io, "Technical options:")
    println(
        io,
        "  --threads <n>             Worker threads to use (default: 1; runtime must provide them)",
    )
    println(io, "  --jobs <n>                Deprecated alias for --threads")
    println(io, "  --quiet                   Suppress informational output")
    println(io, "  --verbose                 Verbose diagnostics to stderr")
    return nothing
end

function _read_model_collection(path::AbstractString, model_type::AbstractString)
    if isdir(path)
        # Directory: load all matching files
        if model_type == "pwm"
            files = filter(f -> endswith(lowercase(f), ".meme"), readdir(path; join=true))
        elseif model_type == "bamm"
            files = filter(f -> endswith(lowercase(f), ".ihbcp"), readdir(path; join=true))
        elseif model_type == "sitega"
            files = filter(f -> endswith(lowercase(f), ".mat"), readdir(path; join=true))
        elseif model_type in ("dimont", "slim")
            files = filter(f -> endswith(lowercase(f), ".xml"), readdir(path; join=true))
        else
            throw(CLIError("unsupported model type for directory: $(model_type)"))
        end
        isempty(files) && throw(CLIError("no $(model_type) files found in $(path)."))
        models = AbstractProfileSource[]
        for f in files
            push!(models, _read_typed_model(f, model_type))
        end
        return models
    else
        # Single file: MEME can contain multiple motifs
        if model_type == "pwm"
            # Read all motifs from MEME
            models = AbstractProfileSource[]
            idx = 0
            while true
                try
                    pfm = read_meme(path; index=idx)
                    pwm = pwm_from_pfm(pfm.frequencies; background=0.25f0, name=pfm.name)
                    push!(models, pwm)
                    idx += 1
                catch e
                    if e isa ModelFormatError && occursin("out of range", e.message)
                        break
                    end
                    rethrow()
                end
            end
            isempty(models) && throw(CLIError("no motifs found in $(path)."))
            return models
        else
            return AbstractProfileSource[_read_typed_model(path, model_type)]
        end
    end
end

function _run_build_null(parsed::CLIParsed)
    if length(parsed.positional) != 1
        _print_build_null_help(stderr)
        throw(CLIError("build-null requires a motif collection path."))
    end
    motifs_path = parsed.positional[1]

    haskey(parsed.options, "model-type") || throw(CLIError("--model-type is required."))
    haskey(parsed.options, "output") || throw(CLIError("--output is required."))
    isdir(motifs_path) || throw(CLIError("motif collection path must be a directory."))

    model_type = parsed.options["model-type"]
    model_type in MODEL_TYPES ||
        throw(CLIError("--model-type must be one of: $(join(MODEL_TYPES, ", "))"))

    output_path = parsed.options["output"]
    seed = _cli_int(parsed, "seed", "127"; minimum=0)
    n_samples = _cli_int(parsed, "num-samples", "2000"; minimum=1)
    num_seq = _cli_int(parsed, "num-sequences", "1000"; minimum=1)
    seq_len = _cli_int(parsed, "seq-length", "200"; minimum=1)
    shuffle = "shuffle" in parsed.flags
    fasta = get(parsed.options, "fasta", nothing)

    metric = get(parsed.options, "metric", "co")
    metric in PROFILE_METRICS ||
        throw(CLIError("--metric must be one of: $(join(PROFILE_METRICS, ", "))"))

    # Read models
    models = _read_model_collection(motifs_path, model_type)

    search_range = tryparse(Int, get(parsed.options, "search-range", "10"))
    window_radius = tryparse(Int, get(parsed.options, "window-radius", "10"))
    realign_window = tryparse(Int, get(parsed.options, "realign-window", "3"))
    min_logfpr = tryparse(Float32, get(parsed.options, "min-logfpr", "0.0"))
    search_range === nothing &&
        throw(CLIError("--search-range must be a non-negative integer."))
    window_radius === nothing &&
        throw(CLIError("--window-radius must be a non-negative integer."))
    realign_window === nothing &&
        throw(CLIError("--realign-window must be a non-negative integer."))
    min_logfpr === nothing && throw(CLIError("--min-logfpr must be a finite number."))
    search_range >= 0 || throw(CLIError("--search-range must be a non-negative integer."))
    window_radius >= 0 || throw(CLIError("--window-radius must be a non-negative integer."))
    realign_window >= 0 ||
        throw(CLIError("--realign-window must be a non-negative integer."))
    isfinite(min_logfpr) || throw(CLIError("--min-logfpr must be a finite number."))

    exec_policy = _execution_policy(parsed)

    sequences = if isnothing(fasta)
        make_random_sequences(num_seq, seq_len; seed=seed)
    else
        first(readsequences(fasta))
    end
    result = build_null(
        models;
        metric=metric,
        n_samples=n_samples,
        shuffle=shuffle,
        seed=seed,
        execution=exec_policy,
        sequences=sequences,
        search_range=search_range,
        window_radius=window_radius,
        realign_window=realign_window,
        min_logfpr=min_logfpr,
    )

    # Save null distribution
    savenull(output_path, result.distribution)

    # Output summary
    summary = Dict{String,Any}(
        "output" => output_path,
        "n_models" => length(models),
        "n_comparisons" => result.total_comparisons,
        "n_null" => result.distribution.n_null,
        "model_type" => result.distribution.model_type,
        "shuffle" => result.distribution.shuffle,
        "seed" => result.distribution.seed,
        "metric" => result.distribution.metric,
        "strategy" => result.distribution.strategy,
    )
    if result.distribution.fit isa GEVFit
        summary["gev_shape"] = result.distribution.fit.shape
        summary["gev_location"] = result.distribution.fit.location
        summary["gev_scale"] = result.distribution.fit.scale
        summary["gev_converged"] = result.distribution.fit.converged
    end
    if "verbose" in parsed.flags
        println(
            stderr,
            "build-null: strategy=$(result.distribution.strategy), metric=$(result.distribution.metric), comparisons=$(result.total_comparisons)",
        )
    end
    "quiet" in parsed.flags || _println_json(summary)
    return 0
end

# ── Command: cache clear ────────────────────────────────────────────────────

function _print_cache_help(io::IO)
    println(io, "Usage: mimosa cache clear [--cache-dir <dir>] [options]")
    println(io, "")
    println(io, "Remove all cached profile artifacts from the specified directory.")
    println(io, "")
    println(io, "Options:")
    println(io, "  --cache-dir <dir>    Cache directory (default: .mimosa-cache)")
    println(io, "  --quiet             Suppress informational output")
    println(io, "  --verbose           Verbose diagnostics to stderr")
    return nothing
end

function _run_cache(parsed::CLIParsed)
    if length(parsed.positional) != 1 || parsed.positional[1] != "clear"
        _print_cache_help(stderr)
        throw(CLIError("cache requires a subcommand: clear"))
    end

    cache_dir = get(parsed.options, "cache-dir", ".mimosa-cache")
    root = abspath(cache_dir)
    (root == dirname(root) || root == homedir()) &&
        throw(CLIError("--cache-dir points to a dangerously broad directory."))
    ispath(root) &&
        (!isdir(root) || islink(root)) &&
        throw(CLIError("--cache-dir must be a real directory, not a file or symlink."))
    cache = Cache(cache_dir)
    removed = clearcache(cache)
    "quiet" in parsed.flags ||
        _println_json(Dict{String,Any}("cache_dir" => cache_dir, "removed" => removed))
    return 0
end

# ── Main entry point ─────────────────────────────────────────────────────────
#
# ArgParse owns CLI syntax, including subcommands, option validation, and help
# output. Runners receive the small adapter struct used by the existing API
# bridge, so scientific validation remains separate from argument parsing.

function _add_common_flags!(settings)
    @add_arg_table! settings begin
        "--quiet"
        action = :store_true
        help = "suppress informational output"
        "--verbose"
        action = :store_true
        help = "enable verbose diagnostics"
    end
end

function _cli_settings()
    settings = ArgParseSettings(;
        prog="mimosa",
        description="Motif comparison and statistical evaluation.",
        version="mimosa $(CLI_VERSION)",
    )
    settings.exit_after_help = false
    settings.exc_handler = (_, err) -> throw(CLIError(err.text))
    @add_arg_table! settings begin
        "--version", "-V"
        action = :show_version
        "profile"
        action = :command
        help = "compare two motif score profiles"
        "build-null"
        action = :command
        help = "build a null distribution from motif comparisons"
        "cache"
        action = :command
        help = "manage the disk cache"
    end

    profile = settings["profile"]
    @add_arg_table! profile begin
        "model1"
        help = "first model or score-profile file"
        required = true
        "model2"
        help = "second model or score-profile file"
        required = true
        "--model1-type"
        required = true
        "--model2-type"
        required = true
        "--metric"
        "--search-range"
        "--window-radius"
        "--realign-window"
        "--min-logfpr"
        "--fasta"
        "--background"
        "--num-sequences"
        "--seq-length"
        "--seed"
        "--background-freq"
        "--threads"
        "--cache-dir"
        "--null-distribution"
        "--effective-number-of-targets"
        "--pvalue"
        action = :store_true
    end
    _add_common_flags!(profile)

    build_null = settings["build-null"]
    @add_arg_table! build_null begin
        "motifs"
        help = "motif collection path"
        required = true
        "--model-type"
        required = true
        "--output"
        required = true
        "--metric"
        "--fasta"
        "--num-sequences"
        "--seq-length"
        "--seed"
        "--num-samples"
        "--search-range"
        "--window-radius"
        "--realign-window"
        "--min-logfpr"
        "--threads"
        "--jobs"
        "--shuffle"
        action = :store_true
    end
    _add_common_flags!(build_null)

    cache = settings["cache"]
    @add_arg_table! cache begin
        "operation"
        help = "cache operation (clear)"
        required = true
        "--cache-dir"
        "--threads"
    end
    _add_common_flags!(cache)
    return settings
end

function _parsed_cli(args::Vector{String})
    parsed_args = ArgParse.parse_args(args, _cli_settings())
    parsed_args === nothing && return nothing
    command = parsed_args["%COMMAND%"]
    command === nothing && throw(CLIError("no command specified."))
    values = parsed_args[command]
    parsed = CLIParsed(command)
    for name in ("model1", "model2", "motifs", "operation")
        haskey(values, name) && push!(parsed.positional, values[name])
    end
    for (name, value) in values
        name in ("model1", "model2", "motifs", "operation") && continue
        if value isa Bool
            value && push!(parsed.flags, name)
        elseif !isnothing(value)
            parsed.options[name] = value
        end
    end
    return parsed
end

"""
    cli_main(args=ARGS)

CLI entry point. Returns an integer exit code (0 on success, 1 on usage error,
2 on runtime error). Prints JSON to stdout on success, diagnostics to stderr.
"""
function cli_main(args::Vector{String}=ARGS)::Int
    try
        parsed = _parsed_cli(args)
        parsed === nothing && return 0
        return _dispatch_runner(parsed.command, parsed)
    catch e
        if e isa CLIError
            println(stderr, "error: $(e.message)")
            return 1
        else
            println(stderr, "error: $(typeof(e).name.name): $(e)")
            return 2
        end
    end
end

function _dispatch_runner(command::AbstractString, parsed::CLIParsed)
    if command == "profile"
        return _run_profile(parsed)
    elseif command == "build-null"
        return _run_build_null(parsed)
    elseif command == "cache"
        return _run_cache(parsed)
    else
        throw(CLIError("unknown command: $(command)"))
    end
end
