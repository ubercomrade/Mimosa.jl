# Strand policies for sequence scanning.

"""
    StrandPolicy

Abstract supertype for strand scanning policies.

Concrete policies:
- [`ForwardOnly`](@ref): scan forward strand only.
- [`ReverseOnly`](@ref): scan reverse complement strand only.
- [`BestStrand`](@ref): scan both strands, return the maximum score per position.
- [`BothStrands`](@ref): scan both strands, return both score tracks.
"""
abstract type StrandPolicy end

"""
    ForwardOnly

Scan the forward strand only.
"""
struct ForwardOnly <: StrandPolicy end

"""
    ReverseOnly

Scan the reverse complement strand only.
"""
struct ReverseOnly <: StrandPolicy end

"""
    BestStrand

Scan both strands and return the maximum score at each position.
"""
struct BestStrand <: StrandPolicy end

"""
    BothStrands

Scan both strands and return both score tracks separately.
"""
struct BothStrands <: StrandPolicy end

"""
    StrandPair{T}

A pair of results for forward and reverse strands.

Fields:
- `forward::T`: forward strand result.
- `reverse::T`: reverse strand result.
"""
struct StrandPair{T}
    forward::T
    reverse::T
end

function Base.show(io::IO, pair::StrandPair)
    return print(io, "StrandPair(forward=$(pair.forward), reverse=$(pair.reverse))")
end
