# Slim (Jstacs) motif type per ADR 0001.
#
# A Slim model stores log-odds scores from a Jstacs GenDisMix classifier
# with a mixture of component and ancestor dependencies. The XML mixture
# parameters (component/ancestor/dependency) are normalized to
# log-probabilities and materialized into a dense 5-ary tensor, then
# flattened to a 2D matrix with axes (context_code, position), identical
# in layout to BaMM and Dimont.
#
# The `span` field is the maximum context dependency: the farthest parent
# position referenced by any motif position. For span=0, the model is
# equivalent to an order-0 BaMM (PWM-like). For span=N, each position may
# depend on up to N preceding positions.
#
# Representation layout:
#   matrix[context_code, position]
#   context_code is a 5-ary integer of length (span + 1)
#   Row indexing: code = b[1] * 5^span + b[2] * 5^(span-1) + ... + b[span+1] * 5^0
#
# Scanning geometry (identical to BaMM/Dimont with order=span):
#   kmer       = span + 1          (bases per scoring term)
#   context    = span              (bases before motif start used for context)
#   window     = motif_len + span  (total sequence window needed)
#   n_terms    = motif_len         (number of scoring terms per window)
#   n_positions = seq_len - window + 1

"""
    Slim{T,M}

Slim motif model from Jstacs (GenDisMix classifier with mixture dependencies).

# Fields
- `name::String`: motif identifier (derived from filename)
- `representation::M`: log-odds matrix with shape `(5^(span+1), motif_length)`,
  indexed by `(context_code, position)`. Row `code` corresponds to the 5-ary
  encoding of `span + 1` consecutive bases.
- `span::Int`: maximum context dependency (0 = independent positions)
- `motif_length::Int`: number of motif positions

The representation is a flattened 2D view of the full `(5, 5, ..., 5,
motif_length)` tensor, materialized from the XML mixture parameters via
log-sum-exp over components and ancestors.
"""
struct Slim{T<:AbstractFloat,M<:AbstractMatrix{T}} <: AbstractMotifModel
    name::String
    representation::M
    span::Int
    motif_length::Int

    function Slim{T,M}(
        name::String, representation::M, span::Int, motif_length::Int
    ) where {T<:AbstractFloat,M<:AbstractMatrix{T}}
        _validate_slim(representation, span, motif_length)
        return new{T,M}(name, representation, span, motif_length)
    end
end

function Slim(
    name::AbstractString,
    representation::AbstractMatrix{T},
    span::Integer,
    motif_length::Integer,
) where {T<:AbstractFloat}
    return Slim{T,typeof(representation)}(
        String(name), representation, Int(span), Int(motif_length)
    )
end

function _validate_slim(representation::AbstractMatrix, span::Int, motif_length::Int)
    if span < 0
        throw(ModelDimensionError("Slim span must be non-negative, got $span."))
    end
    # Guard against exponentiation blow-up: 5^(span+1) rows.
    if span > 10
        throw(
            ModelDimensionError(
                "Slim span must be <= 10 to avoid allocation blow-up, got $span."
            ),
        )
    end
    expected_rows = 5^(span + 1)
    if size(representation, 1) != expected_rows
        throw(
            ModelDimensionError(
                "Slim representation must have $expected_rows rows for span=$span, got $(size(representation, 1)).",
            ),
        )
    end
    if size(representation, 2) != motif_length
        throw(
            ModelDimensionError(
                "Slim representation columns ($(size(representation, 2))) must match motif_length ($motif_length).",
            ),
        )
    end
    if motif_length <= 0
        throw(ModelDimensionError("Slim motif_length must be positive, got $motif_length."))
    end
    if !all(isfinite, representation)
        throw(ModelFormatError("", "Slim representation contains non-finite values."))
    end
    return nothing
end

Base.length(model::Slim) = model.motif_length
Base.eltype(::Type{<:Slim{T}}) where {T} = T
Base.size(model::Slim) = size(model.representation)

function Base.show(io::IO, model::Slim)
    return print(
        io, "Slim(\"$(model.name)\", span=$(model.span), $(size(model.representation)))"
    )
end

function Base.:(==)(a::Slim, b::Slim)
    return a.name == b.name &&
           a.span == b.span &&
           a.motif_length == b.motif_length &&
           a.representation == b.representation
end

function Base.isapprox(a::Slim, b::Slim; kwargs...)
    return a.name == b.name &&
           a.span == b.span &&
           a.motif_length == b.motif_length &&
           isapprox(a.representation, b.representation; kwargs...)
end

"""
    scorebounds(model::Slim)

Return `(min_score, max_score)` theoretical score bounds for a [`Slim`](@ref).

Mirrors Python's `score_bounds_from_representation`: take the per-column min/max
across all context codes and sum across positions.
"""
function scorebounds(model::Slim)
    col_min = vec(minimum(model.representation; dims=1))
    col_max = vec(maximum(model.representation; dims=1))
    return (sum(col_min), sum(col_max))
end

"""
    kmer(model::Slim)

Return the k-mer size (= span + 1) for scanning.
"""
kmer(model::Slim) = model.span + 1

# ── Extensibility API (ADR 0003) ──────────────────────────────────────────────
#
# Slim uses `span` bases preceding the motif site as context. The site
# spans `motif_length` positions; there is no downstream context.

modelname(model::Slim) = model.name
motif_length(model::Slim) = model.motif_length
left_context(model::Slim) = model.span
right_context(::Slim) = 0

"""
    context_length(model::Slim)

Return the context length (= span) for scanning.
"""
context_length(model::Slim) = left_context(model)

"""
    window_size(model::Slim)

Return the total window size needed for scanning (= motif_length + span).
"""
window_size(model::Slim) = model.motif_length + left_context(model)

"""
    scan_width(model::Slim)

Return the number of scanning positions per sequence: `window_size` terms
minus `kmer` plus 1 = `motif_length`.
"""
scan_width(model::Slim) = model.motif_length

"""
    site_start_offset(model::Slim)

Return the offset from scan position to motif start (= `span`): the first
`span` bases of the scan window are context, not motif.
"""
site_start_offset(model::Slim) = left_context(model)
