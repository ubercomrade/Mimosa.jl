# Dimont (Jstacs) motif type per ADR 0001.
#
# A Dimont model stores log-odds scores from a Jstacs Bayesian network
# (MarkovModelDiffSM) with tree-structured context dependencies. The XML
# parameter trees are materialized into a dense 5-ary tensor and flattened
# to a 2D matrix with axes (context_code, position), identical in layout
# to BaMM.
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
# Scanning geometry (identical to BaMM with order=span):
#   kmer       = span + 1          (bases per scoring term)
#   context    = span              (bases before motif start used for context)
#   window     = motif_len + span  (total sequence window needed)
#   n_terms    = motif_len         (number of scoring terms per window)
#   n_positions = seq_len - window + 1

"""
    Dimont{T,M}

Dimont motif model from Jstacs (ThresholdedStrandChIPper with MarkovModelDiffSM).

# Fields
- `name::String`: motif identifier (derived from filename)
- `representation::M`: log-odds matrix with shape `(5^(span+1), motif_length)`,
  indexed by `(context_code, position)`. Row `code` corresponds to the 5-ary
  encoding of `span + 1` consecutive bases.
- `span::Int`: maximum context dependency (0 = independent positions)
- `motif_length::Int`: number of motif positions

The representation is a flattened 2D view of the full `(5, 5, ..., 5,
motif_length)` tensor, materialized from the XML parameter trees.
"""
struct Dimont{T<:AbstractFloat,M<:AbstractMatrix{T}} <: AbstractContextModel{T}
    name::String
    representation::M
    span::Int
    motif_length::Int

    function Dimont{T,M}(
        name::String, representation::M, span::Int, motif_length::Int
    ) where {T<:AbstractFloat,M<:AbstractMatrix{T}}
        _validate_context_model(representation, span, motif_length, "Dimont", "span")
        return new{T,M}(name, representation, span, motif_length)
    end
end

function Dimont(
    name::AbstractString,
    representation::AbstractMatrix{T},
    span::Integer,
    motif_length::Integer,
) where {T<:AbstractFloat}
    return Dimont{T,typeof(representation)}(
        String(name), representation, Int(span), Int(motif_length)
    )
end

function Base.show(io::IO, model::Dimont)
    return print(
        io, "Dimont(\"$(model.name)\", span=$(model.span), $(size(model.representation)))"
    )
end

Base.:(==)(a::Dimont, b::Dimont) = _context_model_equal(a, b)
Base.isapprox(a::Dimont, b::Dimont; kwargs...) = _context_model_isapprox(a, b; kwargs...)

# ── Extensibility API (ADR 0003) ──────────────────────────────────────────────
#
# Dimont uses `span` bases preceding the motif site as context. The site
# spans `motif_length` positions; there is no downstream context.

left_context(model::Dimont) = model.span
