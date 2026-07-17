# Explicit cache for expensive Mimosa computations (scan tracks, comparison
# results, null distributions).
#
# Design principles (PLAN.md Stage 7, REFACTORING.md §15):
#   - No global mutable singleton. The cache is an explicit struct passed
#     around or stored by the caller.
#   - Cache keys are content-based (SHA-256 of serialized parameters), not
#     session-dependent `objectid` or `hash`.
#   - Keys incorporate: algorithm version, schema version, model content
#     fingerprint, config, dtype, and sequence fingerprint.
#   - Atomic writes: temp file + rename. Checksums validated on load.
#   - Corrupted/partial files are treated as cache misses with a diagnostic,
#     never affecting correctness.
#   - Cache is fully disableable (just don't call it).
#   - `import Mimosa` does not create cache directories or touch the filesystem.

using FileWatching: Pidfile
using SHA
using TOML

# ── Schema versions ────────────────────────────────────────────────────────

const CACHE_FORMAT_VERSION = 2
const _CACHE_DATA_NAME = "data.bin"
const _CACHE_META_NAME = "meta.toml"
const _CACHE_KEY_PATTERN = r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$"

# Algorithm version tags: bump when the algorithm changes so stale caches
# are automatically invalidated.
const ALGORITHM_VERSIONS = Dict{String,String}(
    "pwm_scan" => "1",
    "bamm_scan" => "1",
    "sitega_scan" => "1",
    "dimont_scan" => "1",
    "slim_scan" => "1",
    "motif_compare" => "1",
    "profile_compare" => "1",
    "prepared_profile" => "1",
    "null_build" => "1",
)

# This format is intentionally independent of Julia's Serialization stdlib. Cache
# entries can therefore be validated structurally and remain usable across Julia
# sessions that use the same Mimosa cache format.
const PREPARED_PROFILE_CACHE_FORMAT_VERSION = 1
const _PREPARED_PROFILE_CACHE_MAGIC = UInt8[
    0x4d, 0x49, 0x4d, 0x4f, 0x53, 0x41, 0x2d, 0x50, 0x52, 0x45, 0x50, 0x2d, 0x31,
]

# ── Cache type ──────────────────────────────────────────────────────────────

"""
    Cache

An explicit, filesystem-backed cache for expensive Mimosa computations.

Fields:
- `directory::String`: cache root directory.
- `enabled::Bool`: if `false`, all operations are no-ops.

The cache does not create directories at construction time. The directory
is created lazily on the first write. This means `import Mimosa` and
`Cache(dir)` never touch the filesystem.

Create a cache with `Cache(dir)` or `Cache(dir; enabled=false)`. Clear with
[`clearcache`](@ref).
"""
struct Cache
    directory::String
    enabled::Bool
end

"""
    Cache(directory::AbstractString; enabled::Bool=true)

Construct a cache backed by `directory`. The directory is not created until
the first write. Pass `enabled=false` to disable all caching.
"""
function Cache(directory::AbstractString; enabled::Bool=true)
    return Cache(String(directory), enabled)
end

function Base.show(io::IO, cache::Cache)
    return print(io, "Cache(\"$(cache.directory)\", enabled=$(cache.enabled))")
end

# ── Content fingerprinting ──────────────────────────────────────────────────

"""
    content_fingerprint(data::AbstractVector{UInt8})

Return a hex-encoded SHA-256 fingerprint of raw byte data.
"""
function content_fingerprint(data::AbstractVector{UInt8})
    return bytes2hex(SHA.sha256(data))
end

"""
    content_fingerprint(data::AbstractVector{<:Integer})

Return a hex-encoded SHA-256 fingerprint of an integer vector (e.g. offsets).
"""
function content_fingerprint(data::AbstractVector{<:Integer})
    io = IOBuffer()
    write(io, "integer-vector|")
    for x in data
        # Textual canonicalization makes width and signedness explicit and is
        # independent of host byte order.
        write(io, string(typeof(x)), ":", string(x), ";")
    end
    return content_fingerprint(take!(io))
end

"""
    content_fingerprint(s::AbstractString)

Return a hex-encoded SHA-256 fingerprint of a string's UTF-8 bytes.
"""
function content_fingerprint(s::AbstractString)
    return content_fingerprint(Vector{UInt8}(codeunits(s)))
end

"""
    content_fingerprint(arr::AbstractArray{T}) where T<:AbstractFloat

Return a hex-encoded SHA-256 fingerprint of a numeric array. The fingerprint
is based on the element type, dimensions, and raw bytes, making it stable
across Julia sessions.
"""
function content_fingerprint(arr::AbstractArray{T}) where {T<:AbstractFloat}
    io = IOBuffer()
    write(io, string(T))
    write(io, ":")
    for d in size(arr)
        write(io, string(d))
        write(io, ",")
    end
    write(io, ";")
    # Preserve the exact value and width without depending on host endianness.
    for x in arr
        write(io, bitstring(x), ";")
    end
    return content_fingerprint(take!(io))
end

"""
    content_fingerprint(model::AbstractMotifModel)

Return a hex-encoded SHA-256 fingerprint of a motif model's content.
The fingerprint incorporates the model type, name, representation, and
background (for PWM). Context-model fingerprints also include the geometry
contract version so incompatible scan semantics cannot reuse cached artifacts.
"""
function content_fingerprint(model::AbstractMotifModel)
    io = IOBuffer()
    write(io, string(typeof(model)))
    write(io, "|")
    write(io, modelname(model))
    write(io, "|")
    _write_model_fingerprint_body!(io, model)
    return content_fingerprint(take!(io))
end

# Dispatch on concrete built-in types so fingerprint changes are explicit. The
# generic `AbstractMotifModel` method above calls `modelname(model)` and
# `_write_model_fingerprint_body!`, which throws a clear error for unknown
# custom model types.

const _CONTEXT_GEOMETRY_FINGERPRINT = ",geometry=symmetric-v1"

function _write_context_geometry_fingerprint!(io::IO, context::Int)
    context == 0 || write(io, _CONTEXT_GEOMETRY_FINGERPRINT)
    return nothing
end

function _write_model_fingerprint_body!(io::IO, model::PWM)
    write(io, content_fingerprint(model.representation))
    write(io, "|")
    return write(io, join(string.(model.background), ","))
end

function _write_model_fingerprint_body!(io::IO, model::BaMM)
    write(io, content_fingerprint(model.representation))
    write(io, "|")
    write(io, "order=" * string(model.order))
    write(io, ",ml=" * string(model.motif_length))
    return _write_context_geometry_fingerprint!(io, model.order)
end

function _write_model_fingerprint_body!(io::IO, model::SiteGA)
    write(io, content_fingerprint(model.representation))
    write(io, "|")
    return write(io, "ml=" * string(model.motif_length))
end

# Dimont and Slim share the same field layout (`span`, `motif_length`,
# `representation`) and therefore the same fingerprint body. Dispatch
# on each concrete type so a future change to one does not silently
# affect the other.
function _write_model_fingerprint_body!(io::IO, model::Dimont)
    write(io, content_fingerprint(model.representation))
    write(io, "|")
    write(io, "span=" * string(model.span))
    write(io, ",ml=" * string(model.motif_length))
    return _write_context_geometry_fingerprint!(io, model.span)
end

function _write_model_fingerprint_body!(io::IO, model::Slim)
    write(io, content_fingerprint(model.representation))
    write(io, "|")
    write(io, "span=" * string(model.span))
    write(io, ",ml=" * string(model.motif_length))
    return _write_context_geometry_fingerprint!(io, model.span)
end

function _write_model_fingerprint_body!(io::IO, model::AbstractMotifModel)
    return throw(
        ArgumentError(
            "no content fingerprint is defined for $(typeof(model)); " *
            "implement `Mimosa.model_fingerprint(model::MyModel)` for cache/null capability.",
        ),
    )
end

"""
    content_fingerprint(profile::ScoreProfile)

Return a fingerprint of a precomputed profile's name, score data, offsets,
and storage layout. Profiles with the same name but different scores never
share a fingerprint.
"""
function content_fingerprint(profile::ScoreProfile)
    io = IOBuffer()
    write(io, "ScoreProfile|layout=ragged-column-major|dtype=Float32|")
    write(io, modelname(profile))
    write(io, "|data=")
    write(io, content_fingerprint(profile.scores.data))
    write(io, "|offsets=")
    write(io, content_fingerprint(profile.scores.offsets))
    return content_fingerprint(take!(io))
end

"""
    sequence_fingerprint(batch::EncodedSequenceBatch)

Return a hex-encoded SHA-256 fingerprint of an encoded sequence batch,
incorporating the data and offsets.
"""
function sequence_fingerprint(batch::EncodedSequenceBatch)
    io = IOBuffer()
    write(io, "batch:")
    write(io, content_fingerprint(batch.data))
    write(io, "|")
    write(io, content_fingerprint(batch.offsets))
    return content_fingerprint(take!(io))
end

# ── Cache key construction ──────────────────────────────────────────────────

"""
    cache_key(cache::Cache, algorithm::AbstractString, parts::AbstractString...)

Build a stable cache key (hex SHA-256) from an algorithm name, its version,
and content parts. The key is deterministic across Julia sessions.

Returns a 16-character hex string (first 8 bytes of SHA-256).
"""
function cache_key(cache::Cache, algorithm::AbstractString, parts::AbstractString...)
    algo_version = get(ALGORITHM_VERSIONS, algorithm, "0")
    io = IOBuffer()
    write(io, "v=$CACHE_FORMAT_VERSION\n")
    write(io, "algo=$algorithm\n")
    write(io, "algo_ver=$algo_version\n")
    for p in parts
        write(io, p)
        write(io, "\n")
    end
    full_hash = bytes2hex(SHA.sha256(take!(io)))
    return full_hash[1:16]
end

# ── Prepared profile cache ──────────────────────────────────────────────────

"""
    prepared_profile_cache_key(cache, source, sequences=nothing;
                               background=nothing, min_logfpr=0.0f0)

Return the content-addressed key for a normalized, anchor-indexed profile.
The key includes the source, comparison sequences, optional calibration
background, and anchor threshold. `sequences` must be supplied for motif
models; score profiles do not consume sequence batches.
"""
function prepared_profile_cache_key(
    cache::Cache,
    source::AbstractProfileSource,
    sequences::Union{Nothing,EncodedSequenceBatch}=nothing;
    background::Union{Nothing,EncodedSequenceBatch}=nothing,
    min_logfpr::Real=0.0f0,
)
    is_motif = source isa AbstractMotifModel
    is_motif && sequences === nothing && throw(
        ArgumentError("motif prepared-profile cache keys require comparison sequences.")
    )
    is_motif && validate_model(source; capability=:cache)
    source_fingerprint = model_fingerprint(source)
    sequence_part = sequences === nothing ? "sequences=none" : "sequences=$(sequence_fingerprint(sequences))"
    effective_background = is_motif && background === nothing ? sequences : background
    background_part =
        effective_background === nothing ?
        "background=none" : "background=$(sequence_fingerprint(effective_background))"
    threshold = Float32(min_logfpr)
    isfinite(threshold) || throw(ArgumentError("min_logfpr must be finite."))
    return cache_key(
        cache,
        "prepared_profile",
        "source=$source_fingerprint",
        sequence_part,
        background_part,
        "min_logfpr=$(bitstring(threshold))",
    )
end

function _cached_prepared_profile(
    cache::Union{Nothing,Cache},
    source::AbstractProfileSource,
    sequences::Union{Nothing,EncodedSequenceBatch},
    background::Union{Nothing,EncodedSequenceBatch},
    threshold::Float32,
)
    (cache === nothing || !cache.enabled) && return (nothing, nothing)
    key = prepared_profile_cache_key(
        cache, source, sequences; background=background, min_logfpr=threshold
    )
    data = cache_get(cache, key)
    data === nothing && return (key, nothing)
    return (key, _decode_prepared_profile(data))
end

function _store_prepared_profile!(
    cache::Union{Nothing,Cache}, key::Union{Nothing,String}, profile::PreparedProfile
)
    (cache === nothing || !cache.enabled || key === nothing) && return profile
    cache_set(
        cache,
        key,
        _encode_prepared_profile(profile);
        metadata=Dict(
            "algorithm" => "prepared_profile",
            "prepared_profile_format_version" => PREPARED_PROFILE_CACHE_FORMAT_VERSION,
        ),
    )
    return profile
end

function _encode_prepared_profile(profile::PreparedProfile)
    io = IOBuffer()
    write(io, _PREPARED_PROFILE_CACHE_MAGIC)
    _write_cache_u64!(io, PREPARED_PROFILE_CACHE_FORMAT_VERSION)
    _write_cache_string!(io, profile.name)
    _write_cache_f32!(io, profile.min_logfpr)
    _write_cache_ragged!(io, profile.bundle.forward)
    _write_cache_ragged!(io, profile.bundle.reverse)
    _write_cache_anchor_csr!(io, profile.anchors[1])
    _write_cache_anchor_csr!(io, profile.anchors[2])
    return take!(io)
end

function _decode_prepared_profile(data::AbstractVector{UInt8})
    try
        io = IOBuffer(data)
        read(io, UInt8, length(_PREPARED_PROFILE_CACHE_MAGIC)) == _PREPARED_PROFILE_CACHE_MAGIC ||
            return nothing
        _read_cache_u64(io) == PREPARED_PROFILE_CACHE_FORMAT_VERSION || return nothing
        name = _read_cache_string(io)
        threshold = _read_cache_f32(io)
        isfinite(threshold) || return nothing
        forward = _read_cache_ragged(io)
        reverse = _read_cache_ragged(io)
        anchors = (_read_cache_anchor_csr(io), _read_cache_anchor_csr(io))
        eof(io) || return nothing
        return PreparedProfile(name, StrandPair(forward, reverse), anchors, threshold)
    catch
        # A checksum-valid entry can still be from an unknown or malformed
        # prepared-profile format. Treat it as a normal cache miss.
        return nothing
    end
end

function _write_cache_u64!(io::IO, value::Integer)
    value >= 0 || throw(ArgumentError("cache vector length must be non-negative."))
    return write(io, htol(UInt64(value)))
end

function _read_cache_u64(io::IO)
    value = ltoh(read(io, UInt64))
    value <= UInt64(typemax(Int)) || throw(ArgumentError("cached vector is too large."))
    return Int(value)
end

function _write_cache_string!(io::IO, value::String)
    bytes = Vector{UInt8}(codeunits(value))
    _write_cache_u64!(io, length(bytes))
    return write(io, bytes)
end

function _read_cache_string(io::IO)
    count = _read_cache_u64(io)
    count <= bytesavailable(io) || throw(ArgumentError("truncated cached string."))
    return String(read(io, UInt8, count))
end

function _write_cache_f32!(io::IO, value::Float32)
    return write(io, htol(reinterpret(UInt32, value)))
end

function _read_cache_f32(io::IO)
    return reinterpret(Float32, ltoh(read(io, UInt32)))
end

function _write_cache_int_vector!(io::IO, values::AbstractVector{Int})
    _write_cache_u64!(io, length(values))
    for value in values
        write(io, htol(reinterpret(UInt64, Int64(value))))
    end
    return nothing
end

function _read_cache_int_vector(io::IO)
    count = _read_cache_u64(io)
    count <= bytesavailable(io) ÷ sizeof(UInt64) || throw(ArgumentError("truncated cached integer vector."))
    values = Vector{Int}(undef, count)
    for index in eachindex(values)
        value = reinterpret(Int64, ltoh(read(io, UInt64)))
        typemin(Int) <= value <= typemax(Int) || throw(ArgumentError("cached integer is out of range."))
        values[index] = Int(value)
    end
    return values
end

function _write_cache_ragged!(io::IO, ragged::RaggedArray{Float32})
    _write_cache_u64!(io, length(ragged.data))
    for value in ragged.data
        _write_cache_f32!(io, value)
    end
    _write_cache_int_vector!(io, ragged.offsets)
    return nothing
end

function _read_cache_ragged(io::IO)
    count = _read_cache_u64(io)
    count <= bytesavailable(io) ÷ sizeof(UInt32) || throw(ArgumentError("truncated cached score vector."))
    data = Vector{Float32}(undef, count)
    for index in eachindex(data)
        data[index] = _read_cache_f32(io)
        isfinite(data[index]) || throw(ArgumentError("cached scores must be finite."))
    end
    return RaggedArray(data, _read_cache_int_vector(io))
end

function _write_cache_anchor_csr!(io::IO, csr::AnchorCSR)
    _write_cache_int_vector!(io, csr.positions)
    _write_cache_int_vector!(io, csr.offsets)
    return nothing
end

function _read_cache_anchor_csr(io::IO)
    return AnchorCSR(_read_cache_int_vector(io), _read_cache_int_vector(io))
end

# ── Cache file paths ────────────────────────────────────────────────────────

"""
    cache_path(cache::Cache, key::AbstractString)

Return the filesystem path for a cache entry with the given key.
"""
function _validate_cache_key(key::AbstractString)
    value = String(key)
    isempty(value) && throw(ArgumentError("cache key must not be empty."))
    (value == "." || value == "..") &&
        throw(ArgumentError("cache key must be a single path component."))
    occursin('\0', value) && throw(ArgumentError("cache key must not contain NUL."))
    isabspath(value) && throw(ArgumentError("absolute cache keys are not allowed."))
    occursin(r"[/\\]", value) &&
        throw(ArgumentError("cache key must not contain path separators."))
    occursin(r"^[A-Za-z]:", value) &&
        throw(ArgumentError("cache key must not contain a drive prefix."))
    occursin(_CACHE_KEY_PATTERN, value) || throw(
        ArgumentError("cache key must be 1-128 ASCII letters, digits, '.', '_' or '-'.")
    )
    return value
end

function _cache_entry_dir(cache::Cache, key::AbstractString)
    value = _validate_cache_key(key)
    root = abspath(cache.directory)
    path = joinpath(root, value)
    islink(path) && throw(ArgumentError("cache entry path must not be a symlink."))
    root_real = ispath(root) ? realpath(root) : root
    parent_real = ispath(dirname(path)) ? realpath(dirname(path)) : dirname(path)
    relative = relpath(joinpath(parent_real, basename(path)), root_real)
    separator = Sys.iswindows() ? "\\\\" : "/"
    (relative == ".." || startswith(relative, ".." * separator)) &&
        throw(ArgumentError("cache path escapes cache root."))
    return path
end

function cache_path(cache::Cache, key::AbstractString)
    return joinpath(_cache_entry_dir(cache, key), _CACHE_DATA_NAME)
end

function cache_meta_path(cache::Cache, key::AbstractString)
    return joinpath(_cache_entry_dir(cache, key), _CACHE_META_NAME)
end

# ── Cache get/set ───────────────────────────────────────────────────────────

"""
    cache_has(cache::Cache, key::AbstractString)

Check whether a cache entry exists and is valid (not corrupted).
Returns `false` if the cache is disabled or the entry is missing/corrupted.
"""
function cache_has(cache::Cache, key::AbstractString)
    path = cache_path(cache, key)
    meta_path = cache_meta_path(cache, key)
    !cache.enabled && return false
    isfile(path) && isfile(meta_path) || return false
    # Validate checksum
    return _validate_cache_entry(path, meta_path)
end

"""
    cache_get(cache::Cache, key::AbstractString)

Read a cache entry's binary data. Returns `nothing` if the cache is disabled,
the entry is missing, or the checksum does not match.
"""
function cache_get(cache::Cache, key::AbstractString)
    path = cache_path(cache, key)
    meta_path = cache_meta_path(cache, key)
    !cache.enabled && return nothing
    isfile(path) && isfile(meta_path) || return nothing
    _validate_cache_entry(path, meta_path) || return nothing
    return read(path)
end

"""
    cache_get_meta(cache::Cache, key::AbstractString)

Read and parse the metadata (TOML) for a cache entry.
Returns `nothing` if missing or cache disabled.
"""
function cache_get_meta(cache::Cache, key::AbstractString)
    meta_path = cache_meta_path(cache, key)
    !cache.enabled && return nothing
    isfile(meta_path) || return nothing
    try
        return TOML.parsefile(meta_path)
    catch
        return nothing
    end
end

"""
    cache_set(cache::Cache, key::AbstractString, data::AbstractVector{UInt8};
              metadata=Dict{String,Any}())

Write binary data and metadata using a staged sibling directory.
No-op if the cache is disabled. A cross-process lock serializes commits and
clears; the complete entry directory is committed with one atomic rename.
"""
function cache_set(
    cache::Cache,
    key::AbstractString,
    data::AbstractVector{UInt8};
    metadata::Dict=Dict{String,Any}(),
)
    path = cache_path(cache, key)
    !cache.enabled && return nothing

    checksum = bytes2hex(SHA.sha256(data))
    meta = Dict{String,Any}(
        "format_version" => CACHE_FORMAT_VERSION,
        "checksum" => "sha256:$checksum",
        "size" => length(data),
    )
    for (name, value) in metadata
        name in ("format_version", "checksum", "size") && continue
        meta[name] = value
    end

    return _with_cache_lock(cache) do
        parent = abspath(cache.directory)
        stage = mktempdir(parent; prefix=".mimosa-cache-stage-", cleanup=false)
        entry_stage = joinpath(stage, _validate_cache_key(key))
        mkpath(entry_stage)
        try
            staged_path = joinpath(entry_stage, _CACHE_DATA_NAME)
            staged_meta = joinpath(entry_stage, _CACHE_META_NAME)
            write(staged_path, data)
            open(staged_meta, "w") do io
                return TOML.print(io, meta; sorted=true)
            end
            _flush_file(staged_path)
            _flush_file(staged_meta)
            target = _cache_entry_dir(cache, key)
            backup = nothing
            if ispath(target)
                backup = target * ".backup-" * string(rand(UInt))
                mv(target, backup)
            end
            try
                mv(entry_stage, target)
                backup !== nothing && rm(backup; recursive=true, force=true)
            catch
                ispath(target) && rm(target; recursive=true, force=true)
                backup !== nothing && ispath(backup) && mv(backup, target)
                rethrow()
            end
            return path
        finally
            isdir(stage) && rm(stage; recursive=true, force=true)
        end
    end
end

"""
    clearcache(cache::Cache)

Remove all cache entries from the cache directory.
Does not remove the directory itself.
"""
function clearcache(cache::Cache)
    !cache.enabled && return 0
    isdir(cache.directory) || return 0
    return _with_cache_lock(cache) do
        count = 0
        _cleanup_orphan_stages!(cache)
        for key in readdir(cache.directory)
            try
                entry = _cache_entry_dir(cache, key)
                isdir(entry) || continue
                path = cache_path(cache, key)
                meta_path = cache_meta_path(cache, key)
                (isfile(path) || isfile(meta_path)) || continue
                rm(entry; recursive=true, force=true)
                count += 1
            catch error
                error isa ArgumentError || rethrow()
            end
        end
        return count
    end
end

function _cleanup_orphan_stages!(cache::Cache)
    removed = 0
    for entry in readdir(cache.directory; join=true)
        name = basename(entry)
        if startswith(name, ".mimosa-cache-stage-") && isdir(entry)
            rm(entry; recursive=true, force=true)
            removed += 1
        elseif occursin(r"\.backup-[0-9]+$", name) && isdir(entry)
            rm(entry; recursive=true, force=true)
            removed += 1
        end
    end
    return removed
end

"""
    clearcache(cache::Cache, key::AbstractString)

Remove a single cache entry and its metadata.
"""
function clearcache(cache::Cache, key::AbstractString)
    value = _validate_cache_key(key)
    !cache.enabled && return 0
    isdir(cache.directory) || return 0
    return _with_cache_lock(cache) do
        entry = _cache_entry_dir(cache, value)
        isdir(entry) || return 0
        path = cache_path(cache, value)
        meta_path = cache_meta_path(cache, value)
        (isfile(path) || isfile(meta_path)) || return 0
        rm(entry; recursive=true, force=true)
        return 1
    end
end

# ── Internal helpers ───────────────────────────────────────────────────────

const _CACHE_LOCK_NAME = ".mimosa-cache.lock"

function _with_cache_lock(f::F, cache::Cache) where {F}
    root = abspath(cache.directory)
    isdir(root) || mkpath(root)
    # ponytail: one cache-wide lock; use per-key locks only if measured contention warrants it.
    return Pidfile.mkpidlock(f, joinpath(root, _CACHE_LOCK_NAME); stale_age=300)
end

function _validate_cache_entry(path::AbstractString, meta_path::AbstractString)
    try
        meta = TOML.parsefile(meta_path)
        expected = get(meta, "checksum", "")
        startswith(expected, "sha256:") || return false
        expected_hash = expected[8:end]
        actual_hash = bytes2hex(SHA.sha256(read(path)))
        return actual_hash == expected_hash
    catch
        return false
    end
end

function _flush_file(path::AbstractString)
    open(path, "r+") do io
        return flush(io)
    end
    return nothing
end
