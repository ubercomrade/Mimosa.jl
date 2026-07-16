# Model type hierarchy for Mimosa.

include("types.jl")
include("geometry.jl")
include("validation.jl")
include("pwm.jl")
include("bamm.jl")
include("sitega.jl")
include("dimont.jl")
include("slim.jl")

function _validate_context_model(
    representation::AbstractMatrix,
    context::Int,
    motif_length_value::Int,
    model_name::AbstractString,
    context_name::AbstractString,
)
    context >= 0 || throw(
        ModelDimensionError(
            "$model_name $context_name must be non-negative, got $context."
        ),
    )
    context <= 10 || throw(
        ModelDimensionError(
            "$model_name $context_name must be <= 10 to avoid allocation blow-up, got $context.",
        ),
    )
    expected_rows = 5^(context + 1)
    size(representation, 1) == expected_rows || throw(
        ModelDimensionError(
            "$model_name representation must have $expected_rows rows for $context_name=$context, got $(size(representation, 1)).",
        ),
    )
    size(representation, 2) == motif_length_value || throw(
        ModelDimensionError(
            "$model_name representation columns ($(size(representation, 2))) must match motif_length ($motif_length_value).",
        ),
    )
    motif_length_value > 0 || throw(
        ModelDimensionError(
            "$model_name motif_length must be positive, got $motif_length_value."
        ),
    )
    all(isfinite, representation) || throw(
        ModelFormatError("", "$model_name representation contains non-finite values.")
    )
    return nothing
end

Base.length(model::AbstractContextModel) = model.motif_length
Base.eltype(::Type{<:AbstractContextModel{T}}) where {T} = T
Base.size(model::AbstractContextModel) = size(model.representation)

modelname(model::AbstractContextModel) = model.name
motif_length(model::AbstractContextModel) = model.motif_length
kmer(model::AbstractContextModel) = left_context(model) + 1
context_length(model::AbstractContextModel) = left_context(model)
scan_width(model::AbstractContextModel) = model.motif_length

function scorebounds(model::AbstractContextModel)
    col_min = vec(minimum(model.representation; dims=1))
    col_max = vec(maximum(model.representation; dims=1))
    return (sum(col_min), sum(col_max))
end

function _context_model_equal(a::AbstractContextModel, b::AbstractContextModel)
    return a.name == b.name &&
           left_context(a) == left_context(b) &&
           a.motif_length == b.motif_length &&
           a.representation == b.representation
end

function _context_model_isapprox(
    a::AbstractContextModel, b::AbstractContextModel; kwargs...
)
    return a.name == b.name &&
           left_context(a) == left_context(b) &&
           a.motif_length == b.motif_length &&
           isapprox(a.representation, b.representation; kwargs...)
end
