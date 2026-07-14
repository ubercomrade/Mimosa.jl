# Result type for profile comparisons.

"""
    ComparisonResult

Immutable result of a profile comparison.

Fields:
- `query::String`: query model name.
- `target::String`: target model name.
- `score::Float32`: best profile alignment score (higher is better).
- `offset::Int`: offset/shift of the query relative to the target at the best
  alignment. Positive means the query is shifted right.
- `orientation::String`: one of `"++"`, `"+-"`, `"-+"`, `"--"`.
- `metric::String`: canonical profile metric identifier.
- `n_sites::Int`: number of site windows contributing to the score.
"""
struct ComparisonResult
    query::String
    target::String
    score::Float32
    offset::Int
    orientation::String
    metric::String
    n_sites::Int
end

function ComparisonResult(
    query::AbstractString,
    target::AbstractString,
    score,
    offset::Int,
    orientation::AbstractString,
    metric::AbstractString,
    n_sites::Int,
)
    return ComparisonResult(
        String(query),
        String(target),
        Float32(score),
        offset,
        String(orientation),
        String(metric),
        n_sites,
    )
end

function ComparisonResult(
    query::AbstractString,
    target::AbstractString,
    score,
    offset::Int,
    orientation::AbstractString,
    metric::AbstractString,
)
    return ComparisonResult(query, target, score, offset, orientation, metric, 0)
end
