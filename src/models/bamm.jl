# BaMM (Bayesian Markov Model) motif type per ADR 0001.
#
# A BaMM stores log-odds scores for higher-order Markov dependencies between
# nucleotide positions. The representation is a flat matrix with axes
# (context_code, position) where context_code is a 5-ary integer of length
# (order + 1), covering all combinations of A(0), C(1), G(2), T(3), N(4).
#
# The N-state (code containing 4) holds the per-position minimum score,
# matching the Python representation that materializes a 5-ary tensor.

"""
    BaMM{T,M}

Bayesian Markov Model motif with higher-order dependencies.

# Fields
- `name::String`: motif identifier
- `representation::M`: log-odds matrix with shape `(5^(order+1), motif_length)`,
  indexed by `(context_code, position)`. Row `code` corresponds to the 5-ary
  encoding of `order + 1` consecutive bases.
- `order::Int`: Markov order (0 = independent positions, equivalent to PWM)
- `motif_length::Int`: number of motif positions

The representation is a flattened 2D view of the full `(5, 5, ..., 5,
motif_length)` tensor. Code computation:
`code = b[1] * 5^order + b[2] * 5^(order-1) + ... + b[order+1] * 5^0`
where `b[i]` is the 5-ary encoded base.
"""
struct BaMM{T<:AbstractFloat,M<:AbstractMatrix{T}} <: AbstractContextModel{T}
    name::String
    representation::M
    order::Int
    motif_length::Int

    function BaMM{T,M}(
        name::String, representation::M, order::Int, motif_length::Int
    ) where {T<:AbstractFloat,M<:AbstractMatrix{T}}
        _validate_context_model(representation, order, motif_length, "BaMM", "order")
        return new{T,M}(name, representation, order, motif_length)
    end
end

function BaMM(
    name::AbstractString,
    representation::AbstractMatrix{T},
    order::Integer,
    motif_length::Integer,
) where {T<:AbstractFloat}
    return BaMM{T,typeof(representation)}(
        String(name), representation, Int(order), Int(motif_length)
    )
end

function Base.show(io::IO, model::BaMM)
    return print(
        io, "BaMM(\"$(model.name)\", order=$(model.order), $(size(model.representation)))"
    )
end

Base.:(==)(a::BaMM, b::BaMM) = _context_model_equal(a, b)
Base.isapprox(a::BaMM, b::BaMM; kwargs...) = _context_model_isapprox(a, b; kwargs...)

# ── Extensibility API (ADR 0003) ──────────────────────────────────────────────
#
# BaMM uses `order` bases preceding the motif site as context. The site
# itself spans `motif_length` positions, and there is no downstream
# context. `context_length` is kept as an internal alias that delegates
# to `left_context` for the rolling k-mer kernels.

left_context(model::BaMM) = model.order
