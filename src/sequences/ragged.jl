# RaggedArray: flat offset-based ragged storage for variable-length rows.
# Used for both encoded sequences and score profiles.

"""
    RaggedArray{T,V,I}

Flat offset-based ragged array storing variable-length rows in a single
contiguous buffer.

# Invariants (enforced by constructor)
- `offsets[1] == 1` (Julia one-based indexing).
- `offsets` is monotonically non-decreasing.
- `offsets[end] == length(data) + 1`.
- Empty rows are allowed (consecutive equal offsets).
- Number of rows: `length(offsets) - 1`.

See ADR 0002 for the rationale behind this representation.
"""
struct RaggedArray{T,V<:AbstractVector{T},I<:AbstractVector{Int}}
    data::V
    offsets::I

    function RaggedArray{T,V,I}(
        data::V, offsets::I
    ) where {T,V<:AbstractVector{T},I<:AbstractVector{Int}}
        Base.require_one_based_indexing(data, offsets)
        _validate_ragged_offsets(offsets, length(data))
        return new{T,V,I}(data, offsets)
    end
end

function RaggedArray(data::AbstractVector{T}, offsets::AbstractVector{Int}) where {T}
    return RaggedArray{T,typeof(data),typeof(offsets)}(data, offsets)
end

function _validate_ragged_offsets(offsets::AbstractVector{Int}, data_len::Int)
    isempty(offsets) && throw(ArgumentError("offsets must not be empty."))
    offsets[1] != 1 && throw(ArgumentError("offsets[1] must be 1, got $(offsets[1])."))
    for i in 2:length(offsets)
        offsets[i] < offsets[i - 1] &&
            throw(ArgumentError("offsets must be non-decreasing at index $i."))
    end
    offsets[end] != data_len + 1 && throw(
        ArgumentError(
            "offsets[end] must be length(data)+1=$data_len+1, got $(offsets[end])."
        ),
    )
    return nothing
end

"""
    nrows(rag::RaggedArray)

Return the number of rows in a [`RaggedArray`](@ref).
"""
nrows(rag::RaggedArray) = length(rag.offsets) - 1

"""
    rowlength(rag::RaggedArray, i::Int)

Return the length of row `i` in a [`RaggedArray`](@ref).
"""
rowlength(rag::RaggedArray, i::Int) = rag.offsets[i + 1] - rag.offsets[i]

"""
    row(rag::RaggedArray, i::Int)

Return a zero-copy view of row `i` in a [`RaggedArray`](@ref).
"""
function row(rag::RaggedArray, i::Int)
    len = rowlength(rag, i)
    if len == 0
        return view(rag.data, 1:0)
    end
    return @view rag.data[rag.offsets[i]:(rag.offsets[i + 1] - 1)]
end

Base.length(rag::RaggedArray) = nrows(rag)
Base.firstindex(::RaggedArray) = 1
Base.lastindex(rag::RaggedArray) = nrows(rag)
Base.getindex(rag::RaggedArray, i::Int) = row(rag, i)
Base.eltype(::Type{<:RaggedArray{T}}) where {T} = T
Base.IteratorSize(::Type{<:RaggedArray}) = Base.HasLength()
Base.IteratorEltype(::Type{<:RaggedArray}) = Base.EltypeUnknown()

function Base.iterate(rag::RaggedArray, state::Int=1)
    state > nrows(rag) && return nothing
    return (row(rag, state), state + 1)
end

function Base.:(==)(a::RaggedArray, b::RaggedArray)
    return a.offsets == b.offsets && a.data == b.data
end

function Base.isapprox(a::RaggedArray, b::RaggedArray; kwargs...)
    return a.offsets == b.offsets && isapprox(a.data, b.data; kwargs...)
end

function Base.show(io::IO, rag::RaggedArray)
    return print(
        io,
        "RaggedArray{$(eltype(rag))} with $(nrows(rag)) rows, $(length(rag.data)) elements",
    )
end

"""
    build_ragged(rows::AbstractVector{<:AbstractVector})

Build a [`RaggedArray`](@ref) from a vector of row vectors.
"""
function build_ragged(rows::AbstractVector{<:AbstractVector{T}}) where {T}
    n = length(rows)
    offsets = Vector{Int}(undef, n + 1)
    offsets[1] = 1
    for i in 1:n
        offsets[i + 1] = offsets[i] + length(rows[i])
    end
    data = Vector{T}(undef, offsets[end] - 1)
    for i in 1:n
        r = rows[i]
        dest_start = offsets[i]
        for j in eachindex(r)
            data[dest_start + j - 1] = r[j]
        end
    end
    return RaggedArray(data, offsets)
end

"""
    empty_ragged(::Type{T}) where T

Return an empty [`RaggedArray`](@ref) with zero rows.
"""
function empty_ragged(::Type{T}) where {T}
    return RaggedArray(Vector{T}(), [1])
end
