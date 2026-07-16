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
struct Slim{T<:AbstractFloat,M<:AbstractMatrix{T}} <: AbstractContextModel{T}
    name::String
    representation::M
    span::Int
    motif_length::Int

    function Slim{T,M}(
        name::String, representation::M, span::Int, motif_length::Int
    ) where {T<:AbstractFloat,M<:AbstractMatrix{T}}
        _validate_context_model(representation, span, motif_length, "Slim", "span")
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

function Base.show(io::IO, model::Slim)
    return print(
        io, "Slim(\"$(model.name)\", span=$(model.span), $(size(model.representation)))"
    )
end

Base.:(==)(a::Slim, b::Slim) = _context_model_equal(a, b)
Base.isapprox(a::Slim, b::Slim; kwargs...) = _context_model_isapprox(a, b; kwargs...)

# ── Extensibility API (ADR 0003) ──────────────────────────────────────────────
#
# Slim uses `span` bases preceding the motif site as context. The site
# spans `motif_length` positions; there is no downstream context.

left_context(model::Slim) = model.span
