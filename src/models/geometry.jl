# Public geometry contract for motif models (ADR 0003, Extensibility API Plan §3, §5).
#
# This file defines the *public* geometry accessors that every model must
# satisfy. The minimal contract for a custom model that needs no context
# is to subtype `AbstractMotifModel` and define:
#
#   modelname(model)
#   motif_length(model)
#   scan_kernel!(forward, reverse, model, sequence, n_positions)
#
# `left_context` and `right_context` default to zero. Models with context
# override only the nonzero side(s).
#
# Derived quantities (`window_size`, `npositions`, `site_start_offset`)
# are computed from the three primary accessors. Concrete built-in types
# may keep their own overrides as a temporary migration, but the
# identities in ADR 0003 must hold for every subtype.

"""
    modelname(model)::AbstractString

Return a stable, non-empty instance name for `model`.

`modelname` is part of the minimal model contract: comparison,
null-distribution construction, and result serialization use it as the
displayed query/target name. Algorithms must not require a `name`
field; third-party models implement this method instead of exposing a
field.
"""
function modelname end

"""
    left_context(model)::Integer

Return the number of bases preceding the motif site that are needed to
compute one forward-strand score. Default `0` for models without upstream
context.

Forward and reverse scores at the same scan index share the same physical
window and motif interval. A context model therefore declares the downstream
bases needed for reverse scoring through [`right_context`](@ref).
"""
left_context(model::AbstractMotifModel) = 0

"""
    right_context(model)::Integer

Return the number of bases following the motif site that are needed to
compute one score. Default `0` for models without downstream context.
"""
right_context(model::AbstractMotifModel) = 0

"""
    window_size(model)::Int

Return the full window size needed to score one scan position:
`left_context(model) + motif_length(model) + right_context(model)`.

Built-in context models use symmetric physical windows: their left context
supports forward scoring and their right context supports reverse scoring.
Concrete overrides must satisfy the identity above; the default formula is
preferred for new model types.
"""
function window_size(model::AbstractMotifModel)
    left = Int(left_context(model))
    motif = Int(motif_length(model))
    right = Int(right_context(model))
    return Base.Checked.checked_add(Base.Checked.checked_add(left, motif), right)
end

"""
    npositions(model, sequence_length)::Int

Return the number of scan positions in a sequence of length
`sequence_length`: `max(sequence_length - window_size(model) + 1, 0)`.

Returns zero for empty or too-short sequences. Throws `ArgumentError` for a
negative sequence length.
"""
function npositions(model::AbstractMotifModel, sequence_length::Int)
    return npositions(sequence_length, window_size(model))
end

"""
    site_start_offset(model)::Int

Return the offset from the scan position to the motif start:
`left_context(model)`. Extracted sites span
`[scan_position + site_start_offset(model), ..., + motif_length(model) - 1]`.
"""
site_start_offset(model::AbstractMotifModel) = Int(left_context(model))

# ── Fingerprint capability accessor (Extensibility API Plan §5.3) ────────────
#
# `model_fingerprint` is the public capability used by cache and null
# bundles. It is *not* required for plain `compare`. The default method
# delegates to the existing content-based fingerprint via the storage
# layer; built-in fingerprints change when score-affecting contracts change.

"""
    model_fingerprint(model)::String

Return a stable SHA-256 hex string for `model`. Required only when the
model participates in cache keys or null-distribution compatibility
tracking; plain `compare` does not call this method.

Built-in fingerprints are stable within a scoring contract and change when
score-affecting contracts change. A custom model implements this method only if
it needs cache or null capabilities.
"""
function model_fingerprint end
