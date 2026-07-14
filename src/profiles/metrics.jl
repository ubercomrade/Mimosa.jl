# Profile comparison metric types and dispatch helpers.

"""
    AbstractProfileMetric

Abstract supertype for profile comparison metrics.

Concrete metrics:
- [`OverlapCoefficient`](@ref): pooled CO (intersection / min(sum1, sum2)).
- [`OverlapCoefficientRowwise`](@ref): mean per-window CO.
- [`DiceSimilarity`](@ref): pooled Dice (2·intersection / (sum1 + sum2)).
- [`DiceSimilarityRowwise`](@ref): mean per-window Dice.
- [`CosineSimilarityProfile`](@ref): mean per-window cosine similarity.
"""
abstract type AbstractProfileMetric end

"""
    OverlapCoefficient

Pooled overlap coefficient (CO): `intersection / min(sum1, sum2)` over all
selected windows. Higher is better.
"""
struct OverlapCoefficient <: AbstractProfileMetric end

"""
    OverlapCoefficientRowwise

Row-wise overlap coefficient: mean of per-window CO values.
Higher is better. Windows with zero denominator are excluded from the mean.
"""
struct OverlapCoefficientRowwise <: AbstractProfileMetric end

"""
    DiceSimilarity

Pooled Dice similarity: `2·intersection / (sum1 + sum2)` over all
selected windows. Higher is better.
"""
struct DiceSimilarity <: AbstractProfileMetric end

"""
    DiceSimilarityRowwise

Row-wise Dice similarity: mean of per-window Dice values.
"""
struct DiceSimilarityRowwise <: AbstractProfileMetric end

"""
    CosineSimilarityProfile

Row-wise cosine similarity: mean of per-window cosine similarity.
"""
struct CosineSimilarityProfile <: AbstractProfileMetric end

"""Return the stable CLI/storage name for a profile comparison metric."""
metric_name(::OverlapCoefficient) = "co"
metric_name(::OverlapCoefficientRowwise) = "co_rowwise"
metric_name(::DiceSimilarity) = "dice"
metric_name(::DiceSimilarityRowwise) = "dice_rowwise"
metric_name(::CosineSimilarityProfile) = "cosine"

"""
    parse_profile_metric(name)

Convert a profile metric string (`co`, `co_rowwise`, `dice`, `dice_rowwise`,
`cosine`) to a typed metric value.
"""
function parse_profile_metric(name::AbstractString)
    name == "co" && return OverlapCoefficient()
    name == "co_rowwise" && return OverlapCoefficientRowwise()
    name == "dice" && return DiceSimilarity()
    name == "dice_rowwise" && return DiceSimilarityRowwise()
    name == "cosine" && return CosineSimilarityProfile()
    return throw(
        ArgumentError(
            "profile metric must be one of: 'co', 'co_rowwise', 'dice', 'dice_rowwise', 'cosine', got '$name'.",
        ),
    )
end

function _resolve_profile_metric(metric)
    metric isa AbstractProfileMetric && return metric
    metric isa Symbol && return parse_profile_metric(string(metric))
    return parse_profile_metric(metric)
end

"""Return whether a metric pools values over all selected windows."""
is_pooled(::AbstractProfileMetric) = false
is_pooled(::OverlapCoefficient) = true
is_pooled(::DiceSimilarity) = true

"""Return whether a metric uses the Dice denominator."""
is_dice_metric(::AbstractProfileMetric) = false
is_dice_metric(::DiceSimilarity) = true
is_dice_metric(::DiceSimilarityRowwise) = true

const PROFILE_EPS = 1e-6
