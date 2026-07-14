# JSON serialization of comparison results (schema v1 draft).

using Printf

"""Version of the JSON contract used for CLI significance annotations."""
const ANNOTATED_RESULT_SCHEMA_VERSION = 1

"""
    to_dict(result::ComparisonResult)

Return the public dictionary payload for a [`ComparisonResult`](@ref), matching
the Python CLI JSON schema. Significance fields are omitted when not present.
"""
function to_dict(result::ComparisonResult)
    d = Dict{String,Any}(
        "query" => result.query,
        "target" => result.target,
        "score" => Float64(result.score),
        "offset" => result.offset,
        "orientation" => result.orientation,
        "metric" => result.metric,
    )
    if result.n_sites > 0
        d["n_sites"] = result.n_sites
    end
    return d
end

"""
    to_dict(result::AnnotatedResult)

Return the public dictionary payload for an [`AnnotatedResult`](@ref),
including significance fields (`p-value`, `adj.p-value`, `E-value`, `null_id`,
`null_n`, `null_estimator`) when present.
"""
function to_dict(result::AnnotatedResult)
    d = Dict{String,Any}(
        "annotation_schema_version" => ANNOTATED_RESULT_SCHEMA_VERSION,
        "query" => result.query,
        "target" => result.target,
        "score" => Float64(result.score),
        "offset" => result.offset,
        "orientation" => result.orientation,
        "metric" => result.metric,
    )
    if result.n_sites > 0
        d["n_sites"] = result.n_sites
    end
    if result.p_value !== nothing
        d["p-value"] = result.p_value
    end
    if result.adj_p_value !== nothing
        d["adj.p-value"] = result.adj_p_value
    end
    if result.e_value !== nothing
        d["E-value"] = result.e_value
    end
    if result.null_id !== nothing
        d["null_id"] = result.null_id
    end
    if result.null_n !== nothing
        d["null_n"] = result.null_n
    end
    if result.null_estimator !== nothing
        d["null_estimator"] = result.null_estimator
    end
    return d
end

# Minimal JSON value escaping for strings and numbers, producing output
# identical to Python's json.dumps with sort_keys and indent=2.
const _json_escape_chars = Dict(
    '"' => "\\\"", '\\' => "\\\\", '\n' => "\\n", '\r' => "\\r", '\t' => "\\t"
)

function _json_string(s::AbstractString)
    io = IOBuffer()
    write(io, '"')
    for c in s
        if haskey(_json_escape_chars, c)
            write(io, _json_escape_chars[c])
        elseif Int(c) < 0x20
            write(io, @sprintf("\\u%04x", Int(c)))
        else
            write(io, c)
        end
    end
    write(io, '"')
    return String(take!(io))
end

_json_number(x::Real) = repr(x)
_json_number(x::Float64) = _json_float(x)

function _json_float(x::Float64)
    if isinf(x)
        return x > 0 ? "Infinity" : "-Infinity"
    elseif isnan(x)
        return "NaN"
    end
    return repr(x)
end

"""
    to_json(result::Union{ComparisonResult,AnnotatedResult})

Serialize a comparison result to a JSON string matching the Python CLI output.
Keys are written in alphabetical order with 2-space indentation, no trailing
newline, matching `json.dumps(result.to_dict(), indent=2, sort_keys=True)`.
"""
function to_json(result::Union{ComparisonResult,AnnotatedResult})
    d = to_dict(result)
    keys_sorted = sort!(collect(keys(d)))
    io = IOBuffer()
    write(io, "{\n")
    n = length(keys_sorted)
    for (i, k) in enumerate(keys_sorted)
        write(io, "  ")
        write(io, _json_string(k))
        write(io, ": ")
        v = d[k]
        if v isa AbstractString
            write(io, _json_string(v))
        elseif v isa AbstractFloat
            write(io, _json_float(Float64(v)))
        else
            write(io, _json_number(v))
        end
        if i < n
            write(io, ",")
        end
        write(io, "\n")
    end
    write(io, "}")
    return String(take!(io))
end
