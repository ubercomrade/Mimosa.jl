# SiteGA dinucleotide motif type per ADR 0001.
#
# A SiteGA model stores scores for dinucleotide (pair) dependencies between
# adjacent positions. The representation is a flat matrix with axes
# (dinucleotide_code, position) where dinucleotide_code is a 5-ary integer
# of length 2, covering all combinations of A(0), C(1), G(2), T(3), N(4).
#
# Code computation: `code = base1 * 5 + base2` (matches Python's C-order reshape
# of (5, 5, length) → (25, length)).
#
# Scanning geometry (distinct from BaMM):
#   kmer       = 2           (dinucleotide, fixed)
#   context    = 0           (no context before window)
#   n_terms    = motif_length - 1   (one fewer term than positions)
#   window     = motif_length       (window equals motif length)
#   positions  = seq_len - motif_length + 1

"""
    SiteGA{T,M}

Dinucleotide motif model from SiteGA discovery tool.

# Fields
- `name::String`: motif identifier
- `representation::M`: score matrix with shape `(25, motif_length)`, indexed by
  `(dinucleotide_code, position)`. Row `code` corresponds to the 5-ary encoding
  `base1 * 5 + base2` of two consecutive bases.
- `motif_length::Int`: number of motif positions

The representation is a flattened 2D view of the `(5, 5, motif_length)` tensor.
Only the first `motif_length - 1` columns are used in scanning (dinucleotide
terms), but all `motif_length` columns participate in score bounds.
"""
struct SiteGA{T<:AbstractFloat,M<:AbstractMatrix{T}} <: AbstractMotifModel
    name::String
    representation::M
    motif_length::Int

    function SiteGA{T,M}(
        name::String, representation::M, motif_length::Int
    ) where {T<:AbstractFloat,M<:AbstractMatrix{T}}
        _validate_sitega(representation, motif_length)
        return new{T,M}(name, representation, motif_length)
    end
end

function SiteGA(
    name::AbstractString, representation::AbstractMatrix{T}, motif_length::Integer
) where {T<:AbstractFloat}
    return SiteGA{T,typeof(representation)}(String(name), representation, Int(motif_length))
end

function _validate_sitega(representation::AbstractMatrix, motif_length::Int)
    if size(representation, 1) != 25
        throw(
            ModelDimensionError(
                "SiteGA representation must have 25 rows (5×5 dinucleotides), got $(size(representation, 1)).",
            ),
        )
    end
    if size(representation, 2) != motif_length
        throw(
            ModelDimensionError(
                "SiteGA representation columns ($(size(representation, 2))) must match motif_length ($motif_length).",
            ),
        )
    end
    if motif_length <= 0
        throw(
            ModelDimensionError("SiteGA motif_length must be positive, got $motif_length.")
        )
    end
    if !all(isfinite, representation)
        throw(ModelFormatError("", "SiteGA representation contains non-finite values."))
    end
    return nothing
end

Base.length(model::SiteGA) = model.motif_length
Base.eltype(::Type{<:SiteGA{T}}) where {T} = T
Base.size(model::SiteGA) = size(model.representation)

function Base.show(io::IO, model::SiteGA)
    return print(io, "SiteGA(\"$(model.name)\", $(size(model.representation)))")
end

function Base.:(==)(a::SiteGA, b::SiteGA)
    return a.name == b.name &&
           a.motif_length == b.motif_length &&
           a.representation == b.representation
end

function Base.isapprox(a::SiteGA, b::SiteGA; kwargs...)
    return a.name == b.name &&
           a.motif_length == b.motif_length &&
           isapprox(a.representation, b.representation; kwargs...)
end

"""
    scorebounds(model::SiteGA)

Return `(min_score, max_score)` theoretical score bounds for a [`SiteGA`](@ref).

Mirrors Python's `score_bounds_from_representation`: take the per-column min/max
across all dinucleotide codes and sum across positions.
"""
function scorebounds(model::SiteGA)
    col_min = vec(minimum(model.representation; dims=1))
    col_max = vec(maximum(model.representation; dims=1))
    return (sum(col_min), sum(col_max))
end

"""
    kmer(model::SiteGA)

Return the k-mer size (= 2 for dinucleotide SiteGA).
"""
kmer(::SiteGA) = 2

# ── Extensibility API (ADR 0003) ──────────────────────────────────────────────
#
# SiteGA scores adjacent dinucleotides inside the motif window. There is
# no context before or after the site.

modelname(model::SiteGA) = model.name
motif_length(model::SiteGA) = model.motif_length
left_context(::SiteGA) = 0
right_context(::SiteGA) = 0

"""
    context_length(model::SiteGA)

Return the context length (= 0 for SiteGA: no context before the window).
"""
context_length(::SiteGA) = 0

"""
    window_size(model::SiteGA)

Return the total window size needed for scanning (= motif_length).
"""
window_size(model::SiteGA) = model.motif_length

"""
    scan_width(model::SiteGA)

Return the number of scoring terms per window (= motif_length - 1).
"""
scan_width(model::SiteGA) = model.motif_length - 1

"""
    site_start_offset(model::SiteGA)

Return the offset from scan position to motif start (= 0): SiteGA has no
context before the motif window.
"""
site_start_offset(::SiteGA) = 0
