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
struct BaMM{T<:AbstractFloat,M<:AbstractMatrix{T}} <: AbstractMotifModel
    name::String
    representation::M
    order::Int
    motif_length::Int

    function BaMM{T,M}(
        name::String, representation::M, order::Int, motif_length::Int
    ) where {T<:AbstractFloat,M<:AbstractMatrix{T}}
        _validate_bamm(representation, order, motif_length)
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

function _validate_bamm(representation::AbstractMatrix, order::Int, motif_length::Int)
    if order < 0
        throw(ModelDimensionError("BaMM order must be non-negative, got $order."))
    end
    # Guard against exponentiation blow-up: 5^(order+1) rows.
    if order > 10
        throw(
            ModelDimensionError(
                "BaMM order must be <= 10 to avoid allocation blow-up, got $order."
            ),
        )
    end
    expected_rows = 5^(order + 1)
    if size(representation, 1) != expected_rows
        throw(
            ModelDimensionError(
                "BaMM representation must have $expected_rows rows for order=$order, got $(size(representation, 1)).",
            ),
        )
    end
    if size(representation, 2) != motif_length
        throw(
            ModelDimensionError(
                "BaMM representation columns ($(size(representation, 2))) must match motif_length ($motif_length).",
            ),
        )
    end
    if motif_length <= 0
        throw(ModelDimensionError("BaMM motif_length must be positive, got $motif_length."))
    end
    if !all(isfinite, representation)
        throw(ModelFormatError("", "BaMM representation contains non-finite values."))
    end
    return nothing
end

Base.length(model::BaMM) = model.motif_length
Base.eltype(::Type{<:BaMM{T}}) where {T} = T
Base.size(model::BaMM) = size(model.representation)

function Base.show(io::IO, model::BaMM)
    return print(
        io, "BaMM(\"$(model.name)\", order=$(model.order), $(size(model.representation)))"
    )
end

function Base.:(==)(a::BaMM, b::BaMM)
    return a.name == b.name &&
           a.order == b.order &&
           a.motif_length == b.motif_length &&
           a.representation == b.representation
end

function Base.isapprox(a::BaMM, b::BaMM; kwargs...)
    return a.name == b.name &&
           a.order == b.order &&
           a.motif_length == b.motif_length &&
           isapprox(a.representation, b.representation; kwargs...)
end

"""
    scorebounds(model::BaMM)

Return `(min_score, max_score)` theoretical score bounds for a [`BaMM`](@ref).

Mirrors Python's `score_bounds_from_representation`: take the per-column min/max
across all context codes and sum across positions.
"""
function scorebounds(model::BaMM)
    col_min = vec(minimum(model.representation; dims=1))
    col_max = vec(maximum(model.representation; dims=1))
    return (sum(col_min), sum(col_max))
end

"""
    kmer(model::BaMM)

Return the k-mer size (= order + 1) for scanning.
"""
kmer(model::BaMM) = model.order + 1

# ── Extensibility API (ADR 0003) ──────────────────────────────────────────────
#
# BaMM uses `order` bases preceding the motif site as context. The site
# itself spans `motif_length` positions, and there is no downstream
# context. `context_length` is kept as an internal alias that delegates
# to `left_context` for the rolling k-mer kernels.

modelname(model::BaMM) = model.name
motif_length(model::BaMM) = model.motif_length
left_context(model::BaMM) = model.order
right_context(::BaMM) = 0

"""
    context_length(model::BaMM)

Return the context length (= order) for scanning: the number of bases before
the motif start position that contribute to the first term's context.
"""
context_length(model::BaMM) = left_context(model)

"""
    window_size(model::BaMM)

Return the total window size needed for scanning (= motif_length + order).
"""
window_size(model::BaMM) = model.motif_length + left_context(model)

"""
    scan_width(model::BaMM)

Return the number of scanning positions per sequence: `window_size` terms
minus `kmer` plus 1 = `motif_length`.
"""
scan_width(model::BaMM) = model.motif_length

"""
    site_start_offset(model::BaMM)

Return the offset from scan position to motif start (= `order`): the first
`context_length` bases of the scan window are context, not motif.
"""
site_start_offset(model::BaMM) = left_context(model)
