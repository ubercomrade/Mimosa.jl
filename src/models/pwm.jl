# Concrete motif model types per ADR 0001.

const NUCLEOTIDE_CARDINALITY = 4
const PSEUDOCOUNT_PWM::Float32 = 1e-4

"""
    PWM{T,M,B}

Position Weight Matrix: log-odds weights for scanning.

`representation` uses axes `(base, position)` with `base ∈ 1:5` (A, C, G, T, N).
The fifth row holds the N-state score (minimum over concrete bases), matching
the Python representation that materializes a 5-row extended PWM.
`background` is a 4-tuple of nucleotide background frequencies.
"""
struct PWM{T<:AbstractFloat,M<:AbstractMatrix{T},B<:NTuple{4,AbstractFloat}} <:
       AbstractMotifModel
    name::String
    representation::M
    background::B

    function PWM{T,M,B}(
        name::String, representation::M, background::B
    ) where {T<:AbstractFloat,M<:AbstractMatrix{T},B<:NTuple{4,AbstractFloat}}
        _validate_pwm_weights(representation, background)
        return new{T,M,B}(name, representation, background)
    end
end

function PWM(
    name::AbstractString, representation::AbstractMatrix{T}, background::NTuple{4}
) where {T<:AbstractFloat}
    return PWM{T,typeof(representation),typeof(background)}(
        String(name), representation, background
    )
end

function _validate_pwm_weights(weights::AbstractMatrix, background::NTuple{4})
    if size(weights, 1) != 5
        throw(
            ModelDimensionError(
                "PWM weights must have 5 rows (A,C,G,T,N), got $(size(weights, 1))."
            ),
        )
    end
    if size(weights, 2) < 1
        throw(
            ModelDimensionError(
                "PWM motif length must be positive, got $(size(weights, 2))."
            ),
        )
    end
    if !all(isfinite, weights)
        throw(ModelFormatError("", "PWM weights contain non-finite values."))
    end
    # Validate background: finite, non-negative, sum approx 1.
    for i in 1:4
        if !isfinite(Float64(background[i]))
            throw(ModelFormatError("", "PWM background[$i] is not finite."))
        end
        if Float64(background[i]) < 0
            throw(ModelFormatError("", "PWM background[$i] is negative."))
        end
    end
    bg_sum = sum(Float64.(background))
    if !isapprox(bg_sum, 1.0; rtol=1e-4)
        throw(
            ModelFormatError(
                "", "PWM background sum is $bg_sum, expected approximately 1.0."
            ),
        )
    end
    return nothing
end

Base.length(model::PWM) = size(model.representation, 2)

# ── Extensibility API (ADR 0003) ──────────────────────────────────────────────
#
# PWM is a no-context matrix model: `motif_length == length`, both
# contexts are zero, the window equals the motif, and the site starts at
# the scan position.

modelname(model::PWM) = model.name
left_context(::PWM) = 0
right_context(::PWM) = 0

"""
    motif_length(model::AbstractMotifModel)

Return the number of motif positions represented by `model`.
"""
motif_length(model::PWM) = length(model)
window_size(model::PWM) = motif_length(model)

"""
    scorematrix(model::AbstractMotifModel)

Return the matrix used by the scanning kernels. Matrix motifs expose their
frequency or weight matrix; higher-order motifs expose their flattened
context-by-position representation.
"""
function scorematrix end

"""
    scoretype(model::AbstractMotifModel)

Return the element type of [`scorematrix`](@ref) for `model`.
"""
scoretype(model::AbstractMotifModel) = eltype(scorematrix(model))

"""
    site_start_offset(model::PWM)

Return the offset from scan position to motif start (= 0 for PWM/PFM:
no context before the motif window).
"""
site_start_offset(::PWM) = 0

Base.eltype(::Type{<:PWM{T}}) where {T} = T

Base.size(model::PWM) = size(model.representation)

function Base.:(==)(a::PWM, b::PWM)
    return a.name == b.name &&
           a.representation == b.representation &&
           a.background == b.background
end

function Base.isapprox(a::PWM, b::PWM; kwargs...)
    return a.name == b.name &&
           isapprox(a.representation, b.representation; kwargs...) &&
           isapprox(collect(a.background), collect(b.background); kwargs...)
end

"""
    kmer(::PWM)

Return the number of encoded bases per PWM scoring term.
"""
kmer(::PWM) = 1

"""
    context_length(::PWM)

Return the number of bases preceding a PWM motif start used as context.
"""
context_length(::PWM) = 0

"""
    scan_width(model::PWM)

Return the number of PWM scoring terms in one scan window.
"""
scan_width(model::PWM) = motif_length(model)

"""
    npositions(model::PWM, seq_len)

Return the number of PWM scan positions in a sequence of length `seq_len`.
"""
npositions(model::PWM, seq_len::Int) = npositions(seq_len, motif_length(model))

"""
    pfm_to_pwm(pfm; background=0.25)

Convert a Position Frequency Matrix to a log-odds Position Weight Matrix.

The result has 4 rows (base × position), matching the Python `pfm_to_pwm`
which computes `log((pfm + 0.0001) / 0.25)`.
"""
function pfm_to_pwm(
    pfm::AbstractMatrix{T}; background::AbstractFloat=0.25f0
) where {T<:AbstractFloat}
    if size(pfm, 1) != NUCLEOTIDE_CARDINALITY
        throw(ModelDimensionError("PFM must have 4 rows, got $(size(pfm, 1))."))
    end
    pc = T(PSEUDOCOUNT_PWM)
    bg = T(background)
    return @. log((pfm + pc) / bg)
end

"""
    extend_pwm_with_n(weights4)

Extend a 4-row PWM to a 5-row PWM by appending an N-state row equal to the
per-column minimum, matching Python's `pwm_model_from_pfm`.
"""
function extend_pwm_with_n(weights4::AbstractMatrix{T}) where {T<:AbstractFloat}
    if size(weights4, 1) != NUCLEOTIDE_CARDINALITY
        throw(
            ModelDimensionError(
                "PWM weights must have 4 rows to extend, got $(size(weights4, 1))."
            ),
        )
    end
    n_row = vec(minimum(weights4; dims=1))
    return vcat(weights4, reshape(n_row, 1, :))
end

"""
    pwm_from_pfm(pfm; background=0.25, name="")

Build a ready-to-scan [`PWM`](@ref) from a raw position-frequency matrix.

This mirrors Python's `pwm_model_from_pfm`: `pfm_to_pwm` then extend with the
per-column minimum as the N-state row.
"""
function pwm_from_pfm(
    pfm::AbstractMatrix{T}; background::AbstractFloat=0.25f0, name::AbstractString=""
) where {T<:AbstractFloat}
    pwm4 = pfm_to_pwm(pfm; background=background)
    weights = extend_pwm_with_n(pwm4)
    bg = ntuple(_ -> T(background), 4)
    return PWM(name, weights, bg)
end

"""
    reverse_complement(weights)

Return the reverse-complement of a PWM weights matrix.

For the `(base, position)` layout this flips the base rows (A↔T, C↔G) and
reverses the position columns, matching Python's `pwm[::-1, ::-1]`.
"""
function reverse_complement(weights::AbstractMatrix{T}) where {T<:AbstractFloat}
    nrows_val = size(weights, 1)
    if nrows_val == 4
        return reverse(reverse(weights; dims=1); dims=2)
    elseif nrows_val == 5
        rc = similar(weights)
        W = size(weights, 2)
        @inbounds for p in 1:W
            rc[1, p] = weights[4, W - p + 1]
            rc[2, p] = weights[3, W - p + 1]
            rc[3, p] = weights[2, W - p + 1]
            rc[4, p] = weights[1, W - p + 1]
            rc[5, p] = weights[5, W - p + 1]
        end
        return rc
    else
        throw(
            ModelDimensionError("reverse_complement expects 4 or 5 rows, got $nrows_val.")
        )
    end
end

function reverse_complement(model::PWM)
    return PWM(model.name, reverse_complement(model.representation), model.background)
end

"""
    scorebounds(model::PWM)

Return `(min_score, max_score)` theoretical score bounds for a [`PWM`](@ref).

Mirrors Python's `score_bounds_from_representation`: take the per-column min/max
across all rows and sum across positions.
"""
function scorebounds(model::PWM)
    w = model.representation
    col_min = vec(minimum(w; dims=1))
    col_max = vec(maximum(w; dims=1))
    return (sum(col_min), sum(col_max))
end
