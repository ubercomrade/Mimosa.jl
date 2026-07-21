# Unified rolling k-mer scanning kernel shared by PWM and higher-order motifs.

function npositions(seq_len::Int, motif_width::Int)
    seq_len < 0 && throw(ArgumentError("sequence length must be non-negative."))
    motif_width < 1 && throw(ArgumentError("motif width must be positive."))
    return max(seq_len - motif_width + 1, 0)
end

@inline function _require_scannable(model::AbstractMotifModel)
    is_scannable(model) ||
        throw(ArgumentError("$(typeof(model)) is not directly scannable."))
    return nothing
end

function _validate_scan_input(seq::AbstractVector{UInt8}, n_pos::Int, width::Int, dests...)
    Base.require_one_based_indexing(seq)
    for dest in dests
        Base.require_one_based_indexing(dest)
    end
    expected_n_pos = npositions(length(seq), width)
    n_pos == expected_n_pos || throw(
        ArgumentError(
            "n_pos=$n_pos does not match sequence geometry; expected $expected_n_pos for width=$width.",
        ),
    )
    any(code -> code > N_CODE, seq) && throw(ArgumentError("invalid encoded DNA code."))
    any(length(dest) < n_pos for dest in dests) &&
        throw(ArgumentError("destination is too short."))
    return nothing
end

function _scan_dest(data::AbstractVector, offsets::Vector{Int}, row_index::Int)
    start = offsets[row_index]
    stop = offsets[row_index + 1] - 1
    return start > stop ? view(data, 1:0) : @view(data[start:stop])
end

# ── Rolling k-mer code preparation ──────────────────────────────────────────

@inline function _ho_oriented_base(
    seq::AbstractVector{UInt8}, index::Int, reverse_complement::Bool
)
    if !(0 <= index < length(seq))
        return 4
    end
    sequence_index = reverse_complement ? length(seq) - index : index + 1
    base = Int(@inbounds seq[sequence_index])
    return reverse_complement && base != 4 ? 3 - base : base
end

"""
    _ho_kmer_codes(seq, kmer_size, first_start, n_codes; reverse_complement=false)

Build `n_codes` consecutive 5-ary k-mer codes beginning at zero-based
`first_start`. Out-of-range bases are N-coded, which is equivalent to scanning
an N-padded sequence. Each subsequent code is derived by removing its leading
base-5 digit and appending the next base.
"""
function _ho_kmer_codes(
    seq::AbstractVector{UInt8},
    kmer_size::Int,
    first_start::Int,
    n_codes::Int;
    reverse_complement::Bool=false,
)
    kmer_size > 0 || throw(ArgumentError("k-mer size must be positive."))
    n_codes >= 0 || throw(ArgumentError("number of k-mer codes must be non-negative."))
    codes = Vector{Int}(undef, n_codes)
    n_codes == 0 && return codes

    code = 0
    @inbounds for offset in 0:(kmer_size - 1)
        code = 5 * code + _ho_oriented_base(seq, first_start + offset, reverse_complement)
    end
    codes[1] = code

    leading_weight = 5^(kmer_size - 1)
    @inbounds for i in 2:n_codes
        start = first_start + i - 2
        code =
            5 *
            (code - _ho_oriented_base(seq, start, reverse_complement) * leading_weight) +
            _ho_oriented_base(seq, start + kmer_size, reverse_complement)
        codes[i] = code
    end
    return codes
end

function _rolling_kmer_scan_codes!(
    forward::AbstractVector{T},
    reverse::AbstractVector{T},
    score_matrix::AbstractMatrix,
    forward_codes::Vector{Int},
    reverse_codes::Vector{Int},
    n_terms::Int,
    n_pos::Int,
) where {T<:AbstractFloat}
    @inbounds for pos in 1:n_pos
        forward_total = zero(T)
        reverse_total = zero(T)
        for term in 0:(n_terms - 1)
            forward_total += score_matrix[forward_codes[pos + term] + 1, term + 1]
            reverse_total += score_matrix[
                reverse_codes[n_pos - pos + term + 1] + 1, term + 1
            ]
        end
        forward[pos] = forward_total
        reverse[pos] = reverse_total
    end
    return (forward, reverse)
end

function _ho_forward_codes(
    model::AbstractMotifModel, seq::AbstractVector{UInt8}, n_pos::Int
)
    n_pos == 0 && return Int[]
    return _ho_kmer_codes(seq, kmer(model), 0, n_pos + scan_width(model) - 1)
end

function _ho_reverse_codes(
    model::AbstractMotifModel, seq::AbstractVector{UInt8}, n_pos::Int
)
    n_pos == 0 && return Int[]
    return _ho_kmer_codes(
        seq, kmer(model), 0, n_pos + scan_width(model) - 1; reverse_complement=true
    )
end

# ── Generic parallel batch scan helpers ─────────────────────────────────────

function _scan_offsets(batch::EncodedSequenceBatch, model::AbstractMotifModel)
    offsets = Vector{Int}(undef, nsequences(batch) + 1)
    offsets[1] = 1
    @inbounds for i in 1:nsequences(batch)
        offsets[i + 1] = offsets[i] + npositions(model, seqlength(batch, i))
    end
    return offsets
end

function _scan_batch(
    model::AbstractMotifModel,
    batch::EncodedSequenceBatch,
    strands::StrandPolicy,
    execution::Execution,
)
    offsets = _scan_offsets(batch, model)
    data = Vector{Float32}(undef, offsets[end] - 1)
    _parallel_for_weighted(execution, diff(offsets)) do i
        dest = _scan_dest(data, offsets, i)
        return _scan_inplace_model!(strands, dest, model, sequence(batch, i), length(dest))
    end
    return RaggedArray(data, offsets)
end

function _scan_batch_both(
    model::AbstractMotifModel, batch::EncodedSequenceBatch, execution::Execution
)
    offsets = _scan_offsets(batch, model)
    forward = Vector{Float32}(undef, offsets[end] - 1)
    reverse = similar(forward)
    _parallel_for_weighted(execution, diff(offsets)) do i
        forward_dest = _scan_dest(forward, offsets, i)
        reverse_dest = _scan_dest(reverse, offsets, i)
        return scan_both!(
            forward_dest, reverse_dest, model, sequence(batch, i), length(forward_dest)
        )
    end
    return StrandPair(RaggedArray(forward, offsets), RaggedArray(reverse, copy(offsets)))
end

# ── Built-in rolling-k-mer backend ──────────────────────────────────────────
#
# The following traits adapt BaMM, SiteGA, Dimont, and Slim to the generic
# scan kernels shared with PWM.
#
# Each model provides the geometry via the trait functions:
#   kmer(model), context_length(model), window_size(model), scan_width(model)
# and stores its scoring matrix in the `representation` field.
# Score loops use precomputed rolling k-mer codes.

const BuiltinRollingModel = Union{PWM,BaMM,SiteGA,Dimont,Slim}

scorematrix(model::BuiltinRollingModel) = model.representation

function scan_kernel!(
    forward::AbstractVector{Float32},
    reverse::AbstractVector{Float32},
    model::BuiltinRollingModel,
    sequence::AbstractVector{UInt8},
    n_positions::Int,
)
    forward_codes = _ho_forward_codes(model, sequence, n_positions)
    reverse_codes = _ho_reverse_codes(model, sequence, n_positions)
    return _rolling_kmer_scan_codes!(
        forward,
        reverse,
        scorematrix(model),
        forward_codes,
        reverse_codes,
        scan_width(model),
        n_positions,
    )
end

# ── Generic single-sequence scan kernels ──────────────────────────────────
#
# Every model implements `scan_kernel!`. Built-in models delegate to the
# shared rolling-k-mer implementation above; custom models provide their own
# implementation through the same extension point.

"""
    scan_forward!(dest, model::AbstractMotifModel, seq, n_pos)

Fill `dest[1:n_pos]` with forward-strand scores for one sequence.
`n_pos` must equal `npositions(model, length(seq))`.

The shared kernel computes both strands; this method returns its forward track.
"""
function scan_forward!(
    dest::AbstractVector{T},
    model::AbstractMotifModel,
    seq::AbstractVector{UInt8},
    n_pos::Int,
) where {T<:AbstractFloat}
    _require_scannable(model)
    _validate_scan_input(seq, n_pos, window_size(model), dest)
    n_pos == 0 && return dest
    fwd = T === Float32 ? dest : Vector{Float32}(undef, n_pos)
    rev = Vector{Float32}(undef, n_pos)
    _scan_kernel_safe!(fwd, rev, model, seq, n_pos)
    T === Float32 || copyto!(dest, 1, fwd, 1, n_pos)
    return dest
end

"""
    scan_reverse!(dest, model::AbstractMotifModel, seq, n_pos)

Fill `dest[1:n_pos]` with reverse-strand scores for one sequence.
`n_pos` must equal `npositions(model, length(seq))`.
"""
function scan_reverse!(
    dest::AbstractVector{T},
    model::AbstractMotifModel,
    seq::AbstractVector{UInt8},
    n_pos::Int,
) where {T<:AbstractFloat}
    _require_scannable(model)
    _validate_scan_input(seq, n_pos, window_size(model), dest)
    n_pos == 0 && return dest
    fwd = Vector{Float32}(undef, n_pos)
    rev = T === Float32 ? dest : Vector{Float32}(undef, n_pos)
    _scan_kernel_safe!(fwd, rev, model, seq, n_pos)
    T === Float32 || copyto!(dest, 1, rev, 1, n_pos)
    return dest
end

"""
    scan_best_strand!(dest, model::AbstractMotifModel, seq, n_pos)

Fill `dest[1:n_pos]` with the per-position maximum of forward and reverse
scores. `n_pos` must equal `npositions(model, length(seq))`. Ties keep the
forward value.
"""
function scan_best_strand!(
    dest::AbstractVector{T},
    model::AbstractMotifModel,
    seq::AbstractVector{UInt8},
    n_pos::Int,
) where {T<:AbstractFloat}
    return best_hits!(dest, model, seq, n_pos)
end

"""
    best_hits!(dest, model::AbstractMotifModel, seq, n_pos)

Compatibility name for [`scan_best_strand!`](@ref).
"""
function best_hits!(
    dest::AbstractVector{T},
    model::AbstractMotifModel,
    seq::AbstractVector{UInt8},
    n_pos::Int,
) where {T<:AbstractFloat}
    _require_scannable(model)
    _validate_scan_input(seq, n_pos, window_size(model), dest)
    n_pos == 0 && return dest
    fwd = Vector{Float32}(undef, n_pos)
    rev = Vector{Float32}(undef, n_pos)
    _scan_kernel_safe!(fwd, rev, model, seq, n_pos)
    @inbounds for i in 1:n_pos
        dest[i] = rev[i] > fwd[i] ? rev[i] : fwd[i]
    end
    return dest
end

"""
    scan_both!(fwd, rev, model::AbstractMotifModel, seq, n_pos)

Fill `fwd` and `rev` with forward and reverse strand scores respectively.
The destinations may have different floating-point element types. `n_pos` must
equal `npositions(model, length(seq))`.
"""
function scan_both!(
    fwd::AbstractVector{F},
    rev::AbstractVector{R},
    model::AbstractMotifModel,
    seq::AbstractVector{UInt8},
    n_pos::Int,
) where {F<:AbstractFloat,R<:AbstractFloat}
    Base.mightalias(fwd, rev) &&
        throw(ArgumentError("forward and reverse destinations must not alias."))
    _require_scannable(model)
    _validate_scan_input(seq, n_pos, window_size(model), fwd, rev)
    n_pos == 0 && return (fwd, rev)
    fwd32 = F === Float32 ? fwd : Vector{Float32}(undef, n_pos)
    rev32 = R === Float32 ? rev : Vector{Float32}(undef, n_pos)
    _scan_kernel_safe!(fwd32, rev32, model, seq, n_pos)
    F === Float32 || copyto!(fwd, 1, fwd32, 1, n_pos)
    R === Float32 || copyto!(rev, 1, rev32, 1, n_pos)
    return (fwd, rev)
end

# ── Generic single-sequence allocating scan ──────────────────────────────

"""
    scan(model::AbstractMotifModel, seq; strands)

Scan a single encoded sequence with a directly scannable motif model.

Returns:
- `Vector{Float32}` for `ForwardOnly`, `ReverseOnly`, `BestStrand`.
- [`StrandPair{Vector{Float32}}`](@ref) for `BothStrands`.
"""
function scan(
    model::AbstractMotifModel,
    seq::AbstractVector{UInt8};
    strands::StrandPolicy=ForwardOnly(),
)
    validate_model(model; capability=:compare)
    _require_scannable(model)
    n_pos = npositions(model, length(seq))
    return _scan_single_model(strands, model, seq, n_pos)
end

function _scan_single_model(
    ::ForwardOnly, model::AbstractMotifModel, seq::AbstractVector{UInt8}, n_pos::Int
)
    dest = Vector{Float32}(undef, n_pos)
    return scan_forward!(dest, model, seq, n_pos)
end

function _scan_single_model(
    ::ReverseOnly, model::AbstractMotifModel, seq::AbstractVector{UInt8}, n_pos::Int
)
    dest = Vector{Float32}(undef, n_pos)
    return scan_reverse!(dest, model, seq, n_pos)
end

function _scan_single_model(
    ::BestStrand, model::AbstractMotifModel, seq::AbstractVector{UInt8}, n_pos::Int
)
    dest = Vector{Float32}(undef, n_pos)
    return scan_best_strand!(dest, model, seq, n_pos)
end

function _scan_single_model(
    ::BothStrands, model::AbstractMotifModel, seq::AbstractVector{UInt8}, n_pos::Int
)
    fwd = Vector{Float32}(undef, n_pos)
    rev = Vector{Float32}(undef, n_pos)
    scan_both!(fwd, rev, model, seq, n_pos)
    return StrandPair(fwd, rev)
end

# ── Generic single-sequence in-place scan ──────────────────────────────────

"""
    scan!(dest, model::AbstractMotifModel, seq; strands)

Fill the first `npositions(model, length(seq))` elements of `dest`; any
remaining elements are unchanged. `BothStrands()` is unsupported because one
destination cannot hold both tracks; use [`scan_both!`](@ref) or allocating
[`scan`](@ref) instead.

Generic method for all directly scannable motif models.
"""
function scan!(
    dest::AbstractVector{T},
    model::AbstractMotifModel,
    seq::AbstractVector{UInt8};
    strands::StrandPolicy=ForwardOnly(),
) where {T<:AbstractFloat}
    validate_model(model; capability=:compare)
    _require_scannable(model)
    n_pos = npositions(model, length(seq))
    if length(dest) < n_pos
        throw(
            ArgumentError("destination has $(length(dest)) elements, need at least $n_pos.")
        )
    end
    return _scan_inplace_model!(strands, dest, model, seq, n_pos)
end

function _scan_inplace_model!(
    ::ForwardOnly, dest::AbstractVector{T}, model::AbstractMotifModel, seq, n_pos
) where {T<:AbstractFloat}
    return scan_forward!(dest, model, seq, n_pos)
end

function _scan_inplace_model!(
    ::ReverseOnly, dest::AbstractVector{T}, model::AbstractMotifModel, seq, n_pos
) where {T<:AbstractFloat}
    return scan_reverse!(dest, model, seq, n_pos)
end

function _scan_inplace_model!(
    ::BestStrand, dest::AbstractVector{T}, model::AbstractMotifModel, seq, n_pos
) where {T<:AbstractFloat}
    return scan_best_strand!(dest, model, seq, n_pos)
end

function _scan_inplace_model!(
    ::BothStrands, dest::AbstractVector{T}, model::AbstractMotifModel, seq, n_pos
) where {T<:AbstractFloat}
    return throw(
        ArgumentError(
            "scan! with BothStrands is not supported; use scan_both! or scan(model, seq; strands=BothStrands()).",
        ),
    )
end

# ── Generic batch scanning (EncodedSequenceBatch) ─────────────────────────

"""
    scan(model::AbstractMotifModel, batch; strands, execution)

Scan all sequences in a batch with a motif model, returning a
[`RaggedArray{Float32}`](@ref) of scores.

For `BothStrands`, returns a [`StrandPair{RaggedArray{Float32}}`](@ref).

With `Execution(n)` for `n > 1`, sequences are processed in parallel at the
top level. Inner scanning kernels remain serial. `Execution(1)` uses the
sequential fast path.

Generic method for all directly scannable motif models.
"""
function _scan_model_batch(
    model::AbstractMotifModel,
    batch::EncodedSequenceBatch;
    strands::StrandPolicy=ForwardOnly(),
    execution::Execution=Execution(),
)
    _require_scannable(model)
    strands isa BothStrands && return _scan_batch_both(model, batch, execution)
    strands isa Union{ForwardOnly,ReverseOnly,BestStrand} ||
        throw(ArgumentError("unsupported strand policy: $(typeof(strands))"))
    return _scan_batch(model, batch, strands, execution)
end

function scan(
    model::AbstractMotifModel,
    batch::EncodedSequenceBatch;
    strands::StrandPolicy=ForwardOnly(),
    execution::Execution=Execution(),
)
    validate_model(model; capability=:compare)
    return _scan_model_batch(model, batch; strands=strands, execution=execution)
end

# ── Generic scan result lengths ─────────────────────────────────────────

"""
    scan_result_lengths(model::AbstractMotifModel, batch)

Return a `Vector{Int}` with the number of scan positions for each sequence.
Generic method for all directly scannable motif models.
"""
function scan_result_lengths(model::AbstractMotifModel, batch::EncodedSequenceBatch)
    validate_model(model; capability=:compare)
    _require_scannable(model)
    return [npositions(model, seqlength(batch, i)) for i in 1:nsequences(batch)]
end
