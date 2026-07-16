# Shared validation and NPY primitives for portable model/null bundles.
#
# This module is the single source of truth for bundle I/O:
#   - Manifest parsing and validation (TOML schema, format version, checksums)
#   - Path safety (traversal, symlink, absolute path prevention)
#   - NPY header parsing and payload reading/writing
#   - Atomic staged-write protocol
#
# `model_storage.jl` (model bundles) and `null_storage.jl` (null distribution
# bundles) both use these shared helpers. The model/null-specific files contain
# only the code that differs between bundle kinds (model-specific shape
# invariants, GEV parameters, compatibility metadata, etc.).
#
# No Dict{String,Any} domain layer is introduced — each reader accesses the
# TOML manifest directly through the typed validation helpers below.

using SHA
using TOML

const BUNDLE_MANIFEST_NAME = "manifest.toml"
const BUNDLE_DATA_DIR = "data"

# These limits apply only to untrusted on-disk bundles. In-memory scientific
# APIs keep their model-specific limits and are not restricted by these values.
const MAX_BUNDLE_MANIFEST_BYTES = 1_048_576
const MAX_BUNDLE_BLOB_BYTES = 1_073_741_824
const MAX_BUNDLE_ARRAYS = 64
const MAX_BUNDLE_RANK = 8
const MAX_BUNDLE_DIMENSION = 100_000_000
const MAX_BUNDLE_ELEMENTS = 100_000_000
const MAX_BUNDLE_ALLOCATION_BYTES = 1_073_741_824
const MAX_NPY_HEADER_BYTES = 1_048_576

const _NPY_MAGIC = UInt8[0x93, 0x4e, 0x55, 0x4d, 0x50, 0x59]

struct _BundleArraySpec
    name::String
    file::String
    dtype::String
    shape::Vector{Int}
    checksum::String
end

struct _NPYHeader
    dtype::String
    shape::Vector{Int}
    fortran_order::Bool
    payload_bytes::Int
    data_offset::Int
end

function _bundle_error(path::AbstractString, message::AbstractString)
    return ModelFormatError(String(path), String(message))
end

function _rethrow_bundle_error(path::AbstractString, err)
    err isa MimosaError && throw(err)
    return throw(_bundle_error(path, "malformed bundle: $(sprint(showerror, err))."))
end

function _required_manifest_string(
    table::AbstractDict, key::AbstractString, path::AbstractString, context::AbstractString
)
    haskey(table, key) || throw(_bundle_error(path, "$context is missing '$key'."))
    value = table[key]
    value isa AbstractString ||
        throw(_bundle_error(path, "$context '$key' must be a string."))
    isempty(value) && throw(_bundle_error(path, "$context '$key' must not be empty."))
    return String(value)
end

function _required_manifest_table(
    table::AbstractDict, key::AbstractString, path::AbstractString, context::AbstractString
)
    haskey(table, key) || throw(_bundle_error(path, "$context is missing '$key'."))
    value = table[key]
    value isa AbstractDict ||
        throw(_bundle_error(path, "$context '$key' must be a TOML table."))
    return value
end

function _required_manifest_bool(
    table::AbstractDict, key::AbstractString, path::AbstractString, context::AbstractString
)
    haskey(table, key) || throw(_bundle_error(path, "$context is missing '$key'."))
    value = table[key]
    value isa Bool || throw(_bundle_error(path, "$context '$key' must be a boolean."))
    return value
end

function _required_manifest_int(
    table::AbstractDict,
    key::AbstractString,
    path::AbstractString,
    context::AbstractString;
    minimum::Integer=typemin(Int),
    maximum::Integer=typemax(Int),
)
    haskey(table, key) || throw(_bundle_error(path, "$context is missing '$key'."))
    value = table[key]
    value isa Integer && !(value isa Bool) ||
        throw(_bundle_error(path, "$context '$key' must be an integer."))
    value < minimum && throw(_bundle_error(path, "$context '$key' must be >= $minimum."))
    value > maximum && throw(_bundle_error(path, "$context '$key' must be <= $maximum."))
    return Int(value)
end

function _required_manifest_shape(value, path::AbstractString, context::AbstractString)
    value isa AbstractVector ||
        throw(_bundle_error(path, "$context shape must be an array."))
    length(value) <= MAX_BUNDLE_RANK ||
        throw(_bundle_error(path, "$context rank exceeds limit $MAX_BUNDLE_RANK."))
    shape = Vector{Int}(undef, length(value))
    for (i, dim) in enumerate(value)
        dim isa Integer && !(dim isa Bool) ||
            throw(_bundle_error(path, "$context dimension $i must be an integer."))
        dim < 0 && throw(_bundle_error(path, "$context dimension $i is negative."))
        dim > MAX_BUNDLE_DIMENSION && throw(
            _bundle_error(
                path, "$context dimension $i exceeds limit $MAX_BUNDLE_DIMENSION."
            ),
        )
        shape[i] = Int(dim)
    end
    return shape
end

function _required_manifest_floats(
    table::AbstractDict,
    key::AbstractString,
    path::AbstractString,
    context::AbstractString;
    expected_length::Union{Nothing,Int}=nothing,
)
    haskey(table, key) || throw(_bundle_error(path, "$context is missing '$key'."))
    value = table[key]
    value isa AbstractVector ||
        throw(_bundle_error(path, "$context '$key' must be an array."))
    expected_length !== nothing &&
        length(value) != expected_length &&
        throw(_bundle_error(path, "$context '$key' must contain $expected_length values."))
    result = Vector{Float64}(undef, length(value))
    for (i, item) in enumerate(value)
        item isa Real && !(item isa Bool) ||
            throw(_bundle_error(path, "$context '$key' value $i must be numeric."))
        converted = try
            Float64(item)
        catch
            throw(_bundle_error(path, "$context '$key' value $i is not representable."))
        end
        isfinite(converted) ||
            throw(_bundle_error(path, "$context '$key' value $i must be finite."))
        result[i] = converted
    end
    return result
end

function _required_manifest_float(
    table::AbstractDict, key::AbstractString, path::AbstractString, context::AbstractString
)
    haskey(table, key) || throw(_bundle_error(path, "$context is missing '$key'."))
    value = table[key]
    value isa Real && !(value isa Bool) ||
        throw(_bundle_error(path, "$context '$key' must be numeric."))
    converted = try
        Float64(value)
    catch
        throw(_bundle_error(path, "$context '$key' is not representable."))
    end
    isfinite(converted) || throw(_bundle_error(path, "$context '$key' must be finite."))
    return converted
end

function _bundle_dtype_bytes(
    dtype::AbstractString, path::AbstractString, context::AbstractString
)
    dtype in ("<f4", "<f8") || throw(
        _bundle_error(
            path, "$context has unsupported dtype '$dtype'; expected '<f4' or '<f8'."
        ),
    )
    return dtype == "<f4" ? 4 : 8
end

function _bundle_shape_payload_bytes(
    shape::AbstractVector{<:Integer},
    dtype::AbstractString,
    path::AbstractString,
    context::AbstractString,
)
    item_size = _bundle_dtype_bytes(dtype, path, context)
    elements = 1
    for (i, dim) in enumerate(shape)
        dim < 0 && throw(_bundle_error(path, "$context dimension $i is negative."))
        dim > MAX_BUNDLE_DIMENSION &&
            throw(_bundle_error(path, "$context dimension $i exceeds bundle limit."))
        if dim == 0
            elements = 0
        elseif elements != 0
            elements > div(MAX_BUNDLE_ELEMENTS, Int(dim)) && throw(
                _bundle_error(path, "$context exceeds the element allocation budget.")
            )
            elements *= Int(dim)
        end
    end
    elements > div(MAX_BUNDLE_ALLOCATION_BYTES, item_size) &&
        throw(_bundle_error(path, "$context exceeds the byte allocation budget."))
    return elements * item_size
end

function _validate_bundle_checksum(
    checksum::AbstractString, path::AbstractString, context::AbstractString
)
    occursin(r"^sha256:[0-9a-f]{64}$", checksum) || throw(
        _bundle_error(path, "$context checksum must match sha256:<64 lowercase hex>.")
    )
    return nothing
end

function _bundle_within_root(root::AbstractString, candidate::AbstractString)
    root_real = realpath(root)
    candidate_real = realpath(candidate)
    relative = relpath(candidate_real, root_real)
    separator = Sys.iswindows() ? "\\" : "/"
    return relative == "." || (relative != ".." && !startswith(relative, ".." * separator))
end

function _resolve_bundle_path(
    root::AbstractString,
    relative::AbstractString,
    path::AbstractString,
    context::AbstractString;
    require_exists::Bool=false,
)
    isempty(relative) && throw(_bundle_error(path, "$context path must not be empty."))
    occursin('\0', relative) && throw(_bundle_error(path, "$context path contains NUL."))
    isabspath(relative) &&
        throw(_bundle_error(path, "$context path must be relative: '$relative'."))
    startswith(relative, "/") &&
        throw(_bundle_error(path, "$context path must be relative: '$relative'."))
    occursin('\\', relative) &&
        throw(_bundle_error(path, "$context path uses unsupported backslash separators."))
    occursin(r"^[A-Za-z]:", relative) &&
        throw(_bundle_error(path, "$context path uses a platform-specific drive prefix."))

    parts = split(relative, '/'; keepempty=true)
    any(isempty, parts) &&
        throw(_bundle_error(path, "$context path contains an empty component."))
    any(part -> part == "." || part == "..", parts) &&
        throw(_bundle_error(path, "$context path contains traversal components."))

    candidate = joinpath(root, relative)
    if ispath(candidate)
        _bundle_within_root(root, candidate) ||
            throw(_bundle_error(path, "$context path escapes the bundle root."))
    elseif require_exists
        throw(_bundle_error(path, "$context file does not exist: '$relative'."))
    else
        parent = dirname(candidate)
        isdir(parent) && (
            _bundle_within_root(root, parent) ||
            throw(_bundle_error(path, "$context path escapes the bundle root."))
        )
    end
    return candidate
end

function _bundle_file_path(
    root::AbstractString,
    relative::AbstractString,
    path::AbstractString,
    context::AbstractString,
)
    candidate = _resolve_bundle_path(root, relative, path, context; require_exists=true)
    isfile(candidate) ||
        throw(_bundle_error(path, "$context is not a regular file: '$relative'."))
    return candidate
end

function _parse_bundle_array(
    root::AbstractString, arrays::AbstractDict, name::AbstractString, path::AbstractString
)
    haskey(arrays, name) || throw(_bundle_error(path, "bundle is missing array '$name'."))
    raw = arrays[name]
    raw isa AbstractDict ||
        throw(_bundle_error(path, "array '$name' must be a TOML table."))
    file = _required_manifest_string(raw, "file", path, "array '$name'")
    dtype = _required_manifest_string(raw, "dtype", path, "array '$name'")
    shape_value = haskey(raw, "shape") ? raw["shape"] : nothing
    shape_value === nothing &&
        throw(_bundle_error(path, "array '$name' is missing 'shape'."))
    shape = _required_manifest_shape(shape_value, path, "array '$name'")
    checksum = _required_manifest_string(raw, "checksum", path, "array '$name'")
    _bundle_dtype_bytes(dtype, path, "array '$name'")
    _validate_bundle_checksum(checksum, path, "array '$name'")
    _resolve_bundle_path(root, file, path, "array '$name'"; require_exists=false)
    _bundle_shape_payload_bytes(shape, dtype, path, "array '$name'")
    return _BundleArraySpec(String(name), file, dtype, shape, checksum)
end

function _read_bundle_manifest(
    root::AbstractString,
    expected_version::Integer;
    expected_kind::Union{Nothing,AbstractString}=nothing,
    manifest_path::AbstractString=joinpath(root, BUNDLE_MANIFEST_NAME),
)
    isdir(root) || throw(_bundle_error(root, "bundle root is not a directory."))
    manifest_rel = if basename(manifest_path) == BUNDLE_MANIFEST_NAME
        BUNDLE_MANIFEST_NAME
    else
        manifest_path
    end
    manifest_file = _bundle_file_path(root, manifest_rel, root, "manifest")
    size = filesize(manifest_file)
    size <= MAX_BUNDLE_MANIFEST_BYTES || throw(
        _bundle_error(
            root, "manifest exceeds size limit $MAX_BUNDLE_MANIFEST_BYTES bytes."
        ),
    )
    manifest = try
        TOML.parsefile(manifest_file)
    catch err
        throw(_bundle_error(root, "invalid TOML manifest: $(sprint(showerror, err))."))
    end
    manifest isa AbstractDict ||
        throw(_bundle_error(root, "manifest must be a TOML table."))

    format = _required_manifest_string(manifest, "format", root, "manifest")
    format == "mimosa" || throw(_bundle_error(root, "unknown bundle format '$format'."))
    version = _required_manifest_int(
        manifest, "format_version", root, "manifest"; minimum=1, maximum=expected_version
    )
    version == expected_version || throw(
        _bundle_error(
            root,
            "unsupported bundle format version $version; expected $expected_version.",
        ),
    )
    kind = _required_manifest_string(manifest, "kind", root, "manifest")
    expected_kind !== nothing &&
        kind != expected_kind &&
        throw(_bundle_error(root, "expected kind '$expected_kind', got '$kind'."))

    arrays = _required_manifest_table(manifest, "arrays", root, "manifest")
    length(arrays) <= MAX_BUNDLE_ARRAYS ||
        throw(_bundle_error(root, "manifest contains too many arrays."))
    total_payload = 0
    for (name, _) in pairs(arrays)
        name isa AbstractString ||
            throw(_bundle_error(root, "array names must be strings."))
        spec = _parse_bundle_array(root, arrays, String(name), root)
        payload = _bundle_shape_payload_bytes(spec.shape, spec.dtype, root, "array '$name'")
        total_payload <= MAX_BUNDLE_ALLOCATION_BYTES - payload ||
            throw(_bundle_error(root, "bundle exceeds the total allocation budget."))
        total_payload += payload
    end
    return manifest
end

function _validate_bundle_blob_size(path::AbstractString)
    size = filesize(path)
    size <= MAX_BUNDLE_BLOB_BYTES || throw(
        _bundle_error(path, "binary blob exceeds size limit $MAX_BUNDLE_BLOB_BYTES bytes."),
    )
    return size
end

function _validate_bundle_array_checksum(
    root::AbstractString, spec::_BundleArraySpec, bundle_path::AbstractString
)
    file_path = _bundle_file_path(root, spec.file, bundle_path, "array '$(spec.name)'")
    _validate_bundle_blob_size(file_path)
    expected_hash = spec.checksum[8:end]
    actual_hash = _file_sha256(file_path)
    actual_hash == expected_hash ||
        throw(_bundle_error(bundle_path, "checksum mismatch for '$(spec.file)'."))
    return file_path
end

function _read_exact(io::IO, count::Int, path::AbstractString, context::AbstractString)
    count >= 0 || throw(_bundle_error(path, "$context has a negative byte count."))
    bytes = read(io, count)
    length(bytes) == count || throw(
        _bundle_error(
            path, "$context is truncated: expected $count bytes, got $(length(bytes))."
        ),
    )
    return bytes
end

function _read_u16le(bytes::AbstractVector{UInt8})
    length(bytes) == 2 || throw(InvariantError("u16 requires two bytes"))
    return Int(bytes[1]) | (Int(bytes[2]) << 8)
end

function _read_u32le(bytes::AbstractVector{UInt8})
    length(bytes) == 4 || throw(InvariantError("u32 requires four bytes"))
    return Int(bytes[1]) | (Int(bytes[2]) << 8) | (Int(bytes[3]) << 16) |
           (Int(bytes[4]) << 24)
end

function _write_u16le(io::IO, value::Integer)
    0 <= value <= typemax(UInt16) || throw(InvariantError("NPY header is too large."))
    write(io, UInt8(value & 0xff), UInt8((value >> 8) & 0xff))
    return nothing
end

function _npy_skip_space(header::AbstractString, index::Int)
    i = index
    while i <= lastindex(header) && header[i] in (' ', '\t', '\r', '\n')
        i = nextind(header, i)
    end
    return i
end

function _npy_expect(
    header::AbstractString, index::Int, expected::Char, path::AbstractString
)
    i = _npy_skip_space(header, index)
    i <= lastindex(header) && header[i] == expected ||
        throw(_bundle_error(path, "NPY header expected '$expected'."))
    return nextind(header, i)
end

function _npy_quoted(header::AbstractString, index::Int, path::AbstractString)
    i = _npy_skip_space(header, index)
    i <= lastindex(header) && header[i] in ('\'', '"') ||
        throw(_bundle_error(path, "NPY header expected a quoted string."))
    quote_char = header[i]
    i = nextind(header, i)
    value = IOBuffer()
    while i <= lastindex(header)
        char = header[i]
        char == quote_char && return String(take!(value)), nextind(header, i)
        char == '\\' && throw(_bundle_error(path, "NPY header escapes are not supported."))
        write(value, char)
        i = nextind(header, i)
    end
    return throw(_bundle_error(path, "NPY header has an unterminated string."))
end

function _npy_token(header::AbstractString, index::Int, path::AbstractString)
    i = _npy_skip_space(header, index)
    start = i
    while i <= lastindex(header) && !(header[i] in (',', ')', '}', ' ', '\t', '\r', '\n'))
        i = nextind(header, i)
    end
    start < i || throw(_bundle_error(path, "NPY header has an empty token."))
    return String(header[start:prevind(header, i)]), i
end

function _npy_shape(header::AbstractString, index::Int, path::AbstractString)
    i = _npy_expect(header, index, '(', path)
    values = Int[]
    i = _npy_skip_space(header, i)
    if i <= lastindex(header) && header[i] == ')'
        return values, nextind(header, i)
    end
    while true
        token, i = _npy_token(header, i, path)
        value = try
            parse(Int, token)
        catch
            throw(_bundle_error(path, "NPY shape contains a non-integer dimension."))
        end
        value >= 0 || throw(_bundle_error(path, "NPY shape contains a negative dimension."))
        push!(values, value)
        i = _npy_skip_space(header, i)
        i <= lastindex(header) || throw(_bundle_error(path, "NPY shape is unterminated."))
        if header[i] == ')'
            return values, nextind(header, i)
        elseif header[i] == ','
            i = nextind(header, i)
            i = _npy_skip_space(header, i)
            if i <= lastindex(header) && header[i] == ')'
                return values, nextind(header, i)
            end
        else
            throw(_bundle_error(path, "NPY shape expected ',' or ')'."))
        end
    end
end

function _parse_npy_header(header::AbstractString, path::AbstractString)
    all(char -> Int(char) <= 0x7f, header) ||
        throw(_bundle_error(path, "NPY header must contain ASCII metadata."))
    i = _npy_expect(header, firstindex(header), '{', path)
    fields = Dict{String,Any}()
    while true
        i = _npy_skip_space(header, i)
        i <= lastindex(header) || throw(_bundle_error(path, "NPY header is unterminated."))
        header[i] == '}' && (i=nextind(header, i); break)
        key, i = _npy_quoted(header, i, path)
        haskey(fields, key) && throw(_bundle_error(path, "NPY header repeats key '$key'."))
        i = _npy_expect(header, i, ':', path)
        i = _npy_skip_space(header, i)
        i <= lastindex(header) ||
            throw(_bundle_error(path, "NPY header has a missing value."))
        value = if header[i] in ('\'', '"')
            parsed, next_index = _npy_quoted(header, i, path)
            i = next_index
            parsed
        elseif header[i] == '('
            parsed, next_index = _npy_shape(header, i, path)
            i = next_index
            parsed
        else
            token, next_index = _npy_token(header, i, path)
            i = next_index
            if token == "True"
                true
            elseif token == "False"
                false
            else
                throw(_bundle_error(path, "NPY header has unsupported value for '$key'."))
            end
        end
        fields[key] = value
        i = _npy_skip_space(header, i)
        i <= lastindex(header) || throw(_bundle_error(path, "NPY header is unterminated."))
        if header[i] == ','
            i = nextind(header, i)
            continue
        elseif header[i] == '}'
            i = nextind(header, i)
            break
        else
            throw(_bundle_error(path, "NPY header expected ',' or '}'."))
        end
    end
    i = _npy_skip_space(header, i)
    i <= lastindex(header) && throw(_bundle_error(path, "NPY header has trailing data."))

    Set(keys(fields)) == Set(("descr", "fortran_order", "shape")) || throw(
        _bundle_error(
            path, "NPY header must contain exactly descr, fortran_order and shape."
        ),
    )
    fields["descr"] isa AbstractString ||
        throw(_bundle_error(path, "NPY descr must be a string."))
    fields["fortran_order"] isa Bool ||
        throw(_bundle_error(path, "NPY fortran_order must be a boolean."))
    fields["shape"] isa Vector{Int} ||
        throw(_bundle_error(path, "NPY shape must be a tuple of integers."))
    return String(fields["descr"]), fields["shape"], fields["fortran_order"]
end

function _read_npy_header(io::IO, path::AbstractString, file_size::Integer)
    _read_exact(io, 6, path, "NPY magic") == _NPY_MAGIC ||
        throw(_bundle_error(path, "NPY magic is invalid."))
    version = _read_exact(io, 2, path, "NPY version")
    major, minor = version
    major in (UInt8(1), UInt8(2)) && minor == UInt8(0) ||
        throw(_bundle_error(path, "unsupported NPY version $major.$minor."))
    length_bytes = _read_exact(io, major == UInt8(1) ? 2 : 4, path, "NPY header length")
    header_length =
        major == UInt8(1) ? _read_u16le(length_bytes) : _read_u32le(length_bytes)
    header_length <= MAX_NPY_HEADER_BYTES ||
        throw(_bundle_error(path, "NPY header exceeds size limit."))
    header_length <= file_size || throw(_bundle_error(path, "NPY header is truncated."))
    header_bytes = _read_exact(io, header_length, path, "NPY header")
    header = try
        String(header_bytes)
    catch err
        throw(
            _bundle_error(path, "NPY header is not valid text: $(sprint(showerror, err))."),
        )
    end
    dtype, shape, fortran_order = _parse_npy_header(header, path)
    payload_bytes = _bundle_shape_payload_bytes(shape, dtype, path, "NPY array")
    data_offset = Int(position(io))
    data_offset % 64 == 0 || throw(_bundle_error(path, "NPY header alignment is invalid."))
    data_offset + payload_bytes <= file_size ||
        throw(_bundle_error(path, "NPY payload is truncated."))
    return _NPYHeader(dtype, shape, fortran_order, payload_bytes, data_offset)
end

function _decode_npy_payload(raw::Vector{UInt8}, info::_NPYHeader, path::AbstractString)
    T = info.dtype == "<f4" ? Float32 : Float64
    values = reinterpret(T, raw)
    if length(info.shape) == 1
        result = Vector{T}(undef, info.shape[1])
        for i in eachindex(result)
            result[i] = ltoh(values[i])
        end
        return result
    elseif length(info.shape) == 2
        nrows, ncols = info.shape
        result = Matrix{T}(undef, nrows, ncols)
        index = 1
        for row in 1:nrows, column in 1:ncols
            result[row, column] = ltoh(values[index])
            index += 1
        end
        return result
    end
    return throw(
        _bundle_error(
            path, "NPY rank $(length(info.shape)) is not supported by this bundle reader."
        ),
    )
end

function _read_npy_array(
    path::AbstractString;
    expected_dtype::Union{Nothing,AbstractString}=nothing,
    expected_shape::Union{Nothing,AbstractVector{<:Integer}}=nothing,
    expected_rank::Union{Nothing,Int}=nothing,
    expected_fortran_order::Union{Nothing,Bool}=false,
)
    file_size = _validate_bundle_blob_size(path)
    open(path, "r") do io
        info = _read_npy_header(io, path, file_size)
        expected_dtype !== nothing &&
            info.dtype != expected_dtype &&
            throw(
                _bundle_error(
                    path, "NPY dtype '$(info.dtype)' does not match '$expected_dtype'."
                ),
            )
        expected_rank !== nothing &&
            length(info.shape) != expected_rank &&
            throw(_bundle_error(path, "NPY rank does not match the manifest."))
        expected_shape !== nothing &&
            info.shape != Int.(expected_shape) &&
            throw(_bundle_error(path, "NPY shape does not match the manifest."))
        expected_fortran_order !== nothing &&
            info.fortran_order != expected_fortran_order &&
            throw(
                _bundle_error(
                    path, "NPY fortran_order does not match the row-major bundle contract."
                ),
            )
        expected_end = info.data_offset + info.payload_bytes
        file_size == expected_end || throw(
            _bundle_error(
                path,
                "NPY payload length mismatch: expected $(info.payload_bytes) bytes.",
            ),
        )
        raw = _read_exact(io, info.payload_bytes, path, "NPY payload")
        return _decode_npy_payload(raw, info, path)
    end
end

function _read_npy_f64(
    path::AbstractString; expected_shape::Union{Nothing,AbstractVector{<:Integer}}=nothing
)
    return _read_npy_array(
        path;
        expected_dtype="<f8",
        expected_shape=expected_shape,
        expected_rank=1,
        expected_fortran_order=false,
    )
end

function _read_raw_f32_2d(
    path::AbstractString; expected_shape::AbstractVector{<:Integer}, expected_bytes::Integer
)
    expected_bytes > 0 || throw(_bundle_error(path, "raw payload must not be empty."))
    file_size = _validate_bundle_blob_size(path)
    file_size == expected_bytes || throw(
        _bundle_error(
            path,
            "raw payload length mismatch: expected $expected_bytes bytes, got $file_size.",
        ),
    )
    _bundle_shape_payload_bytes(expected_shape, "<f4", path, "raw model array") ==
    expected_bytes || throw(_bundle_error(path, "raw payload length disagrees with shape."))
    nrows, ncols = Int.(expected_shape)
    data = Matrix{Float32}(undef, nrows, ncols)
    open(path, "r") do io
        @inbounds for row in 1:nrows, column in 1:ncols
            data[row, column] = reinterpret(Float32, ltoh(read(io, UInt32)))
        end
    end
    all(isfinite, data) ||
        throw(_bundle_error(path, "raw model array contains non-finite values."))
    return data
end

function _npy_shape_text(shape::AbstractVector{<:Integer})
    length(shape) == 1 && return "($(shape[1]),)"
    return "(" * join(shape, ", ") * ")"
end

function _npy_header(dtype::AbstractString, shape::AbstractVector{<:Integer})
    header_dict = "{'descr': '$dtype', 'fortran_order': False, 'shape': $(_npy_shape_text(shape)), }"
    total = 10 + ncodeunits(header_dict) + 1
    padding = (64 - total % 64) % 64
    header = header_dict * repeat(" ", padding) * "\n"
    ncodeunits(header) <= MAX_NPY_HEADER_BYTES ||
        throw(InvariantError("NPY header exceeds the configured size limit."))
    return header
end

function _write_npy(path::AbstractString, data::Vector{Float64})
    header = _npy_header("<f8", [length(data)])
    open(path, "w") do io
        write(io, _NPY_MAGIC, UInt8(1), UInt8(0))
        _write_u16le(io, ncodeunits(header))
        write(io, header)
        for value in data
            write(io, htol(value))
        end
        return flush(io)
    end
    return nothing
end

function _write_raw_f32_2d(path::AbstractString, data::AbstractMatrix{<:AbstractFloat})
    open(path, "w") do io
        for row in axes(data, 1), column in axes(data, 2)
            value = Float32(data[row, column])
            isfinite(value) ||
                throw(InvariantError("model value cannot be represented as Float32."))
            write(io, htol(reinterpret(UInt32, value)))
        end
        return flush(io)
    end
    return nothing
end

function _write_bundle_manifest(path::AbstractString, manifest::AbstractDict)
    open(path, "w") do io
        TOML.print(io, manifest; sorted=true)
        return flush(io)
    end
    return nothing
end

function _file_sha256(path::AbstractString)
    return open(path, "r") do io
        return bytes2hex(SHA.sha256(io))
    end
end

function _with_bundle_write(path::AbstractString, writer::F) where {F}
    target = abspath(String(path))
    if ispath(target)
        throw(InvariantError("bundle target '$target' already exists."))
    end
    parent = dirname(target)
    mkpath(parent)
    prefix = ".$(basename(target)).mimosa-stage-"
    stage = mktempdir(parent; prefix=prefix, cleanup=false)
    mkpath(joinpath(stage, BUNDLE_DATA_DIR))
    try
        writer(target, stage)
        # A directory rename is atomic only when the destination is absent.
        # Refusing overwrite preserves an existing valid bundle on any failure.
        mv(stage, target)
        return target
    catch err
        err isa MimosaError && throw(err)
        throw(
            InvariantError("failed to write bundle '$target': $(sprint(showerror, err)).")
        )
    finally
        isdir(stage) && try
            rm(stage; recursive=true, force=true)
        catch
            # The stage is deliberately left as an orphan for the documented
            # recovery policy when cleanup itself fails.
        end
    end
end

function _with_bundle_write(writer::F, path::AbstractString) where {F<:Function}
    return _with_bundle_write(path, writer)
end
