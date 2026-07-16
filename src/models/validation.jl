# Minimal scan extension point and interface validation
# (Extensibility API Plan §4, §9).
#
# A custom model that needs no context can subtype `AbstractMotifModel`,
# implement `modelname`, `motif_length`, and `scan_kernel!`, and
# automatically participate in scan, prepare_profile, compare, selectsites,
# and reconstruct_pfm without touching Mimosa.jl source.
#
# This file provides:
#   - the safe public wrapper `scan_kernel!` boundary (used by the
#     generic scan fallback),
#   - `validate_model(model; capability=...)` for runtime interface checks,
#   - `ModelInterfaceError` is defined in `errors.jl`.

const SUPPORTED_CAPABILITIES = (:compare, :sites, :cache)

"""
    scan_kernel!(forward, reverse, model, sequence, n_positions)

Fill `forward[1:n_positions]` and `reverse[1:n_positions]` with the
forward and reverse strand scores for `model` against one encoded
sequence, and return `(forward, reverse)`.

This is the *minimal scan extension point* for custom models. The
function is called **after** the public scanning boundary has validated:

- DNA codes and sequence geometry,
- the computed `n_positions`,
- output buffer type, length, and absence of aliasing between
  `forward` and `reverse`.

The kernel must not mutate `model` or `sequence`, must not resize the
buffers, and must write exactly `n_positions` elements to each buffer.

Canonical scanning value type is `Float32`.

Specialized `scan_forward!`, `scan_reverse!`, `best_hits!`, or
`scan_both!` methods may be added as performance overrides, but custom
models must still implement this kernel as their portable scan
capability.
"""
function scan_kernel! end

# ── Specialized-method detection (compile-time capability check) ─────────────
#
# `validate_model` uses these helpers to confirm that a model has a scan
# path. The generic `scan_*!` methods for `AbstractMotifModel` use
# `scan_kernel!`, so detecting "has a scan path" reduces to checking whether
# the model defines the kernel for its concrete type.

function _has_scan_kernel(model)
    return hasmethod(
        scan_kernel!,
        (
            AbstractVector{Float32},
            AbstractVector{Float32},
            typeof(model),
            AbstractVector{UInt8},
            Int,
        ),
    )
end

# ── Safe wrapper that uses the user-provided pair kernel ────────────────────
#
# This wrapper calls model code after the public scanning boundary has
# validated inputs. It enforces the kernel return-value contract and is
# used by `scan_both!`/`scan_forward!`/`scan_reverse!`/`best_hits!`
# defined in `scanning/n_order_scan.jl`.

function _scan_kernel_safe!(
    forward::AbstractVector{Float32},
    reverse::AbstractVector{Float32},
    model::AbstractMotifModel,
    seq::AbstractVector{UInt8},
    n_pos::Int,
)
    result = scan_kernel!(forward, reverse, model, seq, n_pos)
    valid_result =
        result isa Tuple &&
        length(result) == 2 &&
        result[1] === forward &&
        result[2] === reverse
    valid_result || throw(InvariantError("scan_kernel! must return (forward, reverse)."))
    return (forward, reverse)
end

# ── Interface validation ────────────────────────────────────────────────────

function _interface_error(capability::Symbol, model, message::AbstractString)
    return ModelInterfaceError(capability, string(typeof(model)), String(message))
end

function _has_method(f, model, args...)
    return hasmethod(f, (typeof(model), args...))
end

"""
    validate_model(model; capability=:compare) -> model

Validate that `model` satisfies the public extension interface for the
requested `capability`. Returns `model` on success, throws
[`ModelInterfaceError`](@ref) on failure.

## Capabilities

| Capability       | Required methods |
|------------------|------------------|
| `:compare`       | `AbstractMotifModel` subtype, `modelname`, `motif_length`, `scan_kernel!` |
| `:sites`         | `:compare` plus valid geometry: positive `motif_length`, non-negative contexts, `window_size >= motif_length` |
| `:cache`         | `:compare` plus `model_fingerprint` |

Symbols are accepted on this user-facing boundary; internally the
capability value is only used for diagnostics, never for hot-path
dispatch.
"""
function validate_model(model; capability::Symbol=:compare)
    capability in SUPPORTED_CAPABILITIES || throw(
        ArgumentError(
            "capability must be one of $(join(string.(SUPPORTED_CAPABILITIES), ", ")), got :$capability",
        ),
    )

    model isa AbstractMotifModel ||
        throw(_interface_error(capability, model, "model must subtype AbstractMotifModel."))

    # :compare requirements.
    _has_method(modelname, model) || throw(
        _interface_error(
            capability,
            model,
            "missing required method `modelname(model)::AbstractString`.",
        ),
    )
    _has_method(motif_length, model) || throw(
        _interface_error(
            capability, model, "missing required method `motif_length(model)::Integer`."
        ),
    )

    # Geometry validation (shared by :compare and :sites).
    ml = motif_length(model)
    ml_int = _interface_int(ml, capability, model, "motif_length")
    ml_int > 0 || throw(
        _interface_error(capability, model, "motif_length must be positive, got $ml.")
    )
    lc = _interface_int(left_context(model), capability, model, "left_context")
    rc = _interface_int(right_context(model), capability, model, "right_context")
    lc >= 0 || throw(
        _interface_error(capability, model, "left_context must be non-negative, got $lc."),
    )
    rc >= 0 || throw(
        _interface_error(capability, model, "right_context must be non-negative, got $rc."),
    )
    expected_ws = _checked_geometry_sum(lc, ml_int, rc, capability, model)
    actual_ws = _interface_int(window_size(model), capability, model, "window_size")
    actual_ws == expected_ws || throw(
        _interface_error(
            capability,
            model,
            "window_size must equal left_context + motif_length + right_context " *
            "($expected_ws), got $actual_ws.",
        ),
    )

    # Every model declares its scan capability through `scan_kernel!`.
    # Generic `scan_*!` methods are wrappers and do not themselves
    # constitute a scan path.
    _has_scan_kernel(model) || throw(
        _interface_error(
            capability,
            model,
            "missing scan capability: implement `scan_kernel!(forward, reverse, model, seq, n_pos)` " *
            "for the concrete model type.",
        ),
    )

    # name must be a non-empty string.
    name = modelname(model)
    name isa AbstractString ||
        throw(_interface_error(capability, model, "modelname must return a String."))
    isempty(name) &&
        throw(_interface_error(capability, model, "modelname must not be empty."))

    if capability === :sites
        offset = _interface_int(
            site_start_offset(model), capability, model, "site_start_offset"
        )
        offset == lc || throw(
            _interface_error(
                capability,
                model,
                "site_start_offset must equal left_context ($lc), got $offset.",
            ),
        )
        ml_int <= actual_ws - offset || throw(
            _interface_error(capability, model, "motif site lies outside the scan window."),
        )
    end

    if capability === :cache
        fallback = which(model_fingerprint, (AbstractProfileSource,))
        selected = which(model_fingerprint, (typeof(model),))
        selected !== fallback || throw(
            _interface_error(
                capability,
                model,
                "missing `model_fingerprint(model)::String`; required for cache/null compatibility tracking.",
            ),
        )
        fp = model_fingerprint(model)
        fp isa AbstractString || throw(
            _interface_error(capability, model, "model_fingerprint must return a String."),
        )
        isempty(fp) && throw(
            _interface_error(
                capability, model, "model_fingerprint must not return an empty string."
            ),
        )
        occursin(r"^[0-9a-fA-F]{64}$", fp) || throw(
            _interface_error(
                capability,
                model,
                "model_fingerprint must be a 64-character SHA-256 hex string.",
            ),
        )
    end

    return model
end

function _interface_int(value, capability::Symbol, model, accessor::AbstractString)
    value isa Integer ||
        throw(_interface_error(capability, model, "$accessor must return an Integer."))
    try
        converted = Int(value)
        converted == value || throw(InexactError(:Int, Int, value))
        return converted
    catch err
        err isa InexactError || err isa OverflowError || rethrow()
        throw(
            _interface_error(
                capability, model, "$accessor must convert to Int without loss, got $value."
            ),
        )
    end
end

function _checked_geometry_sum(left::Int, motif::Int, right::Int, capability::Symbol, model)
    try
        return Base.Checked.checked_add(Base.Checked.checked_add(left, motif), right)
    catch err
        err isa OverflowError || rethrow()
        throw(
            _interface_error(
                capability,
                model,
                "window_size overflow: $left + $motif + $right does not fit in Int.",
            ),
        )
    end
end
