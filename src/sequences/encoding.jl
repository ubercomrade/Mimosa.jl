# Encoded sequence batch: flat UInt8 buffer with offsets for ragged DNA sequences.
# Per ADR 0002: A=0x00, C=0x01, G=0x02, T=0x03, N/ambiguous/padding=0x04.

"""
    N_CODE

The `UInt8` code used for N (any base) and padding in encoded sequences.
Value is `0x04`. Valid codes are `0x00`–`0x04` (A, C, G, T, N/padding).
"""
const N_CODE = 0x04

# 256-entry lookup table: maps ASCII byte to 5-ary nucleotide code.
# Indexed by byte+1 (Julia 1-based), so _ENCODE_TABLE[byte + 1] gives the code.
const _ENCODE_TABLE = fill(N_CODE, 256)
_ENCODE_TABLE[Int('A') + 1] = 0x00
_ENCODE_TABLE[Int('C') + 1] = 0x01
_ENCODE_TABLE[Int('G') + 1] = 0x02
_ENCODE_TABLE[Int('T') + 1] = 0x03
_ENCODE_TABLE[Int('a') + 1] = 0x00
_ENCODE_TABLE[Int('c') + 1] = 0x01
_ENCODE_TABLE[Int('g') + 1] = 0x02
_ENCODE_TABLE[Int('t') + 1] = 0x03

"""
    encode_base(byte::UInt8)

Encode one ASCII byte to a 5-ary nucleotide code.
A=0, C=1, G=2, T=3, all others (including N and IUPAC) = 4.
"""
encode_base(byte::UInt8) = _ENCODE_TABLE[byte + 1]

"""
    encode_sequence(str::AbstractString)

Encode a DNA string to a `Vector{UInt8}` using 5-ary encoding.
"""
function encode_sequence(str::AbstractString)
    n = ncodeunits(str)
    result = Vector{UInt8}(undef, n)
    # Invariant: _ENCODE_TABLE maps every possible byte (0..255) to a value in
    # 0..N_CODE, so all outputs are valid by construction.  @inbounds is safe
    # because result has exactly n elements and i ranges over 1..n.
    @inbounds for i in 1:n
        result[i] = _ENCODE_TABLE[codeunit(str, i) + 1]
    end
    return result
end

"""
    _validate_encoded_data(data)

Check that every byte in `data` is a valid 5-ary nucleotide code (`0..N_CODE`).
Raises `InvariantError` if any byte is out of range.

This is called at every public construction boundary so that `@inbounds`
scanning kernels can safely use `Int(seq[i]) + 1` as a row index into
representation matrices without runtime bounds checks.
"""
function _validate_encoded_data(data::AbstractVector{UInt8})
    for i in eachindex(data)
        if data[i] > N_CODE
            throw(
                InvariantError(
                    "invalid encoded base 0x$(string(data[i], base=16)) at index $i; " *
                    "valid codes are 0x00..0x04 (A,C,G,T,N).",
                ),
            )
        end
    end
    return nothing
end

"""
    EncodedSequenceBatch{V,I}

A batch of DNA sequences stored as a flat `UInt8` buffer with offsets.

# Invariants
- `offsets[1] == 1`.
- `offsets` is monotonically non-decreasing.
- `offsets[end] == length(data) + 1`.
- Encoding: A=0x00, C=0x01, G=0x02, T=0x03, N/ambiguous=0x04.
- Every byte in `data` satisfies `0 <= byte <= N_CODE`.

All public constructors validate these invariants.  Internal callers that
guarantee validity by construction may use `_unsafe_encoded_batch` to
skip the scan.

See ADR 0002 for the design rationale.
"""
struct EncodedSequenceBatch{V<:AbstractVector{UInt8},I<:AbstractVector{Int}}
    data::V
    offsets::I

    function EncodedSequenceBatch{V,I}(
        data::V, offsets::I
    ) where {V<:AbstractVector{UInt8},I<:AbstractVector{Int}}
        Base.require_one_based_indexing(data, offsets)
        _validate_ragged_offsets(offsets, length(data))
        _validate_encoded_data(data)
        return new{V,I}(data, offsets)
    end

    # Internal unsafe constructor — skips code validation for hot paths
    # where the caller guarantees all bytes are in `0:N_CODE` by construction.
    function EncodedSequenceBatch{V,I}(
        data::V, offsets::I, ::Val{:unsafe}
    ) where {V<:AbstractVector{UInt8},I<:AbstractVector{Int}}
        Base.require_one_based_indexing(data, offsets)
        _validate_ragged_offsets(offsets, length(data))
        return new{V,I}(data, offsets)
    end
end

function EncodedSequenceBatch(data::AbstractVector{UInt8}, offsets::AbstractVector{Int})
    return EncodedSequenceBatch{typeof(data),typeof(offsets)}(data, offsets)
end

# Internal unsafe constructor — skips code validation for hot paths where the
# caller guarantees all bytes are in `0:N_CODE` by construction.
# Used by `make_random_sequences` and `EncodedSequenceBatch(rows::...)`.
function _unsafe_encoded_batch(data::AbstractVector{UInt8}, offsets::AbstractVector{Int})
    _validate_ragged_offsets(offsets, length(data))
    return EncodedSequenceBatch{typeof(data),typeof(offsets)}(data, offsets, Val{:unsafe}())
end

"""
    nsequences(batch::EncodedSequenceBatch)

Return the number of sequences in a batch.
"""
nsequences(batch::EncodedSequenceBatch) = length(batch.offsets) - 1

"""
    seqlength(batch::EncodedSequenceBatch, i::Int)

Return the length of sequence `i` in a batch.
"""
seqlength(batch::EncodedSequenceBatch, i::Int) = batch.offsets[i + 1] - batch.offsets[i]

"""
    sequence(batch::EncodedSequenceBatch, i::Int)

Return a zero-copy view of sequence `i` in a batch.
"""
function sequence(batch::EncodedSequenceBatch, i::Int)
    len = seqlength(batch, i)
    if len == 0
        # Return a safe empty view at valid indices
        return view(batch.data, 1:0)
    end
    return @view batch.data[batch.offsets[i]:(batch.offsets[i + 1] - 1)]
end

Base.length(batch::EncodedSequenceBatch) = nsequences(batch)
Base.firstindex(::EncodedSequenceBatch) = 1
Base.lastindex(batch::EncodedSequenceBatch) = nsequences(batch)
Base.getindex(batch::EncodedSequenceBatch, i::Int) = sequence(batch, i)
Base.IteratorSize(::Type{<:EncodedSequenceBatch}) = Base.HasLength()
Base.IteratorEltype(::Type{<:EncodedSequenceBatch}) = Base.EltypeUnknown()

function Base.iterate(batch::EncodedSequenceBatch, state::Int=1)
    state > nsequences(batch) && return nothing
    return (sequence(batch, state), state + 1)
end

function Base.:(==)(a::EncodedSequenceBatch, b::EncodedSequenceBatch)
    return a.offsets == b.offsets && a.data == b.data
end

function Base.show(io::IO, batch::EncodedSequenceBatch)
    return print(
        io,
        "EncodedSequenceBatch($(nsequences(batch)) sequences, $(length(batch.data)) bytes)",
    )
end

"""
    empty_sequence_batch()

Return an empty [`EncodedSequenceBatch`](@ref) with zero sequences.
"""
empty_sequence_batch() = EncodedSequenceBatch(UInt8[], [1])

# ── Random sequence generation ──────────────────────────────────────────────

# Lookup from a uniform Float64 in [0, 1) to a nucleotide code.
# Bases are assigned equal-width intervals: A=[0,0.25), C=[0.25,0.5),
# G=[0.5,0.75), T=[0.75,1.0).
const _BASE_LOOKUP = (0x00, 0x01, 0x02, 0x03)

"""
    make_random_sequences(rng::AbstractRNG, n::Int, len::Int)

Generate `n` random DNA sequences of length `len` each using `rng`. Bases are
drawn uniformly from A, C, G, T.

Returns an [`EncodedSequenceBatch`](@ref). Reproducible within Julia but
not bit-compatible with Python's `np.random.default_rng` (different RNG
algorithm). This is acceptable for CLI fallback sequences; users should
provide explicit FASTA for scientific reproducibility across languages.
"""
function make_random_sequences(rng::AbstractRNG, n::Int, len::Int)
    n < 0 && throw(ArgumentError("n must be non-negative, got $n."))
    len < 0 && throw(ArgumentError("len must be non-negative, got $len."))
    total = Base.Checked.checked_mul(n, len)
    flat_data = Vector{UInt8}(undef, total)
    offsets = Vector{Int}(undef, n + 1)
    offsets[1] = 1
    for i in 1:n
        start = offsets[i]
        for j in 1:len
            idx = floor(Int, rand(rng) * 4.0) + 1
            flat_data[start + j - 1] = _BASE_LOOKUP[idx]
        end
        offsets[i + 1] = start + len
    end
    # All codes are drawn from _BASE_LOOKUP (0..3), so skip validation.
    return _unsafe_encoded_batch(flat_data, offsets)
end

function make_random_sequences(n::Int, len::Int; seed::Integer=127)
    return make_random_sequences(Random.MersenneTwister(seed), n, len)
end

"""
    EncodedSequenceBatch(rows::AbstractVector{<:AbstractVector{UInt8}})

Build an [`EncodedSequenceBatch`](@ref) from a vector of encoded sequence vectors.
All bytes must satisfy `0 <= byte <= N_CODE`; otherwise an `InvariantError`
is raised.
"""
function EncodedSequenceBatch(rows::AbstractVector{<:AbstractVector{UInt8}})
    n = length(rows)
    offsets = Vector{Int}(undef, n + 1)
    offsets[1] = 1
    for i in 1:n
        offsets[i + 1] = offsets[i] + length(rows[i])
    end
    data = Vector{UInt8}(undef, offsets[end] - 1)
    for i in 1:n
        r = rows[i]
        dest_start = offsets[i]
        for j in eachindex(r)
            data[dest_start + j - 1] = r[j]
        end
    end
    # Validate codes at this public boundary.
    return EncodedSequenceBatch(data, offsets)
end

"""
    reverse_complement(seq::AbstractVector{UInt8})

Return the reverse complement of an encoded DNA sequence.

For 5-ary encoding: complement of A(0)↔T(3), C(1)↔G(2), N(4) stays N(4).
The result is a new vector; the input is not modified.
"""
function reverse_complement(seq::AbstractVector{UInt8})
    n = length(seq)
    result = Vector{UInt8}(undef, n)
    # Invariant: seq codes are valid 0..N_CODE (guaranteed by caller or
    # EncodedSequenceBatch construction).  @inbounds is safe because result
    # has exactly n elements and src indices 1..n are valid.
    @inbounds for i in 1:n
        b = seq[n - i + 1]
        result[i] = b == N_CODE ? N_CODE : 0x03 - b
    end
    return result
end

"""
    reverse_complement!(dest::AbstractVector{UInt8}, src::AbstractVector{UInt8})

Write the reverse complement of `src` into `dest`. Both must have the same length.

# Aliasing
`dest` and `src` must not be the same array.  If they alias, an `ArgumentError`
is raised because the in-place reverse complement would corrupt data during
the copy.
"""
function reverse_complement!(dest::AbstractVector{UInt8}, src::AbstractVector{UInt8})
    n = length(src)
    length(dest) >= n || throw(ArgumentError("dest must be at least as long as src."))
    # Copy before writing whenever views overlap; this also makes identical
    # arrays and partially overlapping views safe and deterministic.
    source = Base.mightalias(dest, src) ? copy(src) : src
    # Invariant: src codes are valid 0..N_CODE (guaranteed by EncodedSequenceBatch
    # construction).  complement(b) = N_CODE if b==N_CODE else 0x03-b, which is
    # also in 0..N_CODE.  @inbounds is safe because dest has >= n elements and
    # src indices 1..n are valid.
    @inbounds for i in 1:n
        b = source[n - i + 1]
        dest[i] = b == N_CODE ? N_CODE : 0x03 - b
    end
    return dest
end

"""
    reverse_complement(batch::EncodedSequenceBatch)

Return a new [`EncodedSequenceBatch`](@ref) where every sequence is reverse-complemented.
"""
function reverse_complement(batch::EncodedSequenceBatch)
    n = nsequences(batch)
    data = Vector{UInt8}(undef, length(batch.data))
    offsets = copy(batch.offsets)
    for i in 1:n
        reverse_complement!(
            @view(data[offsets[i]:(offsets[i + 1] - 1)]), sequence(batch, i)
        )
    end
    return _unsafe_encoded_batch(data, offsets)
end

# Padded conversion helpers for compatibility with oracle fixtures and
# kernel scratch buffers.

"""
    to_padded(batch::EncodedSequenceBatch; padding::UInt8=N_CODE)

Return `(matrix, lengths)` where `matrix` is a dense padded `Matrix{UInt8}`
with `padding` filling unused columns, and `lengths` is the per-sequence length vector.

This is NOT the canonical representation; it is provided for compatibility
testing and kernel scratch buffers.
"""
function to_padded(batch::EncodedSequenceBatch; padding::UInt8=N_CODE)
    n = nsequences(batch)
    n == 0 && return (Matrix{UInt8}(undef, 0, 0), Int[])
    max_len = maximum(i -> seqlength(batch, i), 1:n)
    matrix = fill(padding, n, max_len)
    lengths = Vector{Int}(undef, n)
    for i in 1:n
        len = seqlength(batch, i)
        lengths[i] = len
        for j in 1:len
            matrix[i, j] = batch.data[batch.offsets[i] + j - 1]
        end
    end
    return (matrix, lengths)
end

"""
    from_padded(values::AbstractMatrix{UInt8}, lengths::AbstractVector{Int}; padding::UInt8=N_CODE)

Build an [`EncodedSequenceBatch`](@ref) from a padded dense matrix and a lengths vector.
Only the first `lengths[i]` columns of each row are used.
"""
function from_padded(
    values::AbstractMatrix{UInt8}, lengths::AbstractVector{Int}; padding::UInt8=N_CODE
)
    # Validate that padding is a valid code.
    padding > N_CODE && throw(
        InvariantError("padding must be in 0..N_CODE, got 0x$(string(padding, base=16))."),
    )
    n = length(lengths)
    n == size(values, 1) ||
        throw(ArgumentError("lengths has $n rows, but values has $(size(values, 1))."))
    total = 0
    for i in 1:n
        len = lengths[i]
        len < 0 && throw(ArgumentError("lengths must be non-negative, got $len at row $i."))
        len > size(values, 2) && throw(
            ArgumentError("length $len exceeds matrix width $(size(values, 2)) at row $i."),
        )
        total = Base.Checked.checked_add(total, len)
    end
    data = Vector{UInt8}(undef, total)
    offsets = Vector{Int}(undef, n + 1)
    offsets[1] = 1
    for i in 1:n
        len = lengths[i]
        start = offsets[i]
        @inbounds for j in 1:len
            data[start + j - 1] = values[i, j]
        end
        offsets[i + 1] = start + len
    end
    return EncodedSequenceBatch(data, offsets)
end
