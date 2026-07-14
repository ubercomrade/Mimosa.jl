# Precomputed score profiles: an input to profile comparison, not a motif model.

"""
    ScoreProfile

Precomputed per-position score profiles read from FASTA-like numerical files.
Both strands resolve to the same scores.

Fields:
- `name::String`: profile name (derived from filename).
- `scores::RaggedArray{Float32}`: one row per sequence, variable-length.
"""
struct ScoreProfile <: AbstractProfileSource
    name::String
    scores::RaggedArray{Float32}
end

Base.length(profile::ScoreProfile) = nrows(profile.scores)

"""
    modelname(profile::ScoreProfile)

Return the profile name. ScoreProfile is an `AbstractProfileSource`, not an
`AbstractMotifModel`, but `modelname` is the public accessor for all profile
sources.
"""
modelname(profile::ScoreProfile) = profile.name

function Base.show(io::IO, profile::ScoreProfile)
    return print(io, "ScoreProfile(\"$(profile.name)\", $(nrows(profile.scores)) rows)")
end

function Base.:(==)(a::ScoreProfile, b::ScoreProfile)
    return a.name == b.name && a.scores == b.scores
end

"""
    scorebounds(profile::ScoreProfile)

Return `(min_score, max_score)` from the precomputed score values.
"""
function scorebounds(profile::ScoreProfile)
    if isempty(profile.scores.data)
        return (0.0f0, 0.0f0)
    end
    return (minimum(profile.scores.data), maximum(profile.scores.data))
end

"""
    profile_bundle(profile::ScoreProfile)

Return a `StrandPair{RaggedArray{Float32}}` where both strands are the
same precomputed scores.
"""
function profile_bundle(profile::ScoreProfile)
    return StrandPair(profile.scores, profile.scores)
end
