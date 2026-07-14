# BaMM (.ihbcp) parser: reads Bayesian Markov Model files and converts to
# log-odds representation against a uniform background.
#
# File format:
#   - Blocks separated by blank lines, one block per motif position
#   - Each block has (max_order + 1) lines, one per order k = 0, 1, ..., max_order
#   - Line k has 4^(k+1) conditional probability values
#   - Comments starting with '#' are ignored
#
# The parser converts probabilities to log-odds against a uniform background
# (0.25^(k+1)), handles early positions with truncated order by broadcasting,
# and creates a 5-ary tensor with N-state = per-position minimum.

const MAX_BAMM_POSITIONS = 10_000
const MAX_BAMM_ORDER = 10
const BAMM_EPSILON::Float32 = 1e-10

"""
    read_bamm(path; order=nothing)

Read a BaMM motif from an `.ihbcp` file and return a [`BaMM`](@ref).

If `order` is not specified, the file's maximum order is used. If `order`
exceeds the file's maximum order, the file's maximum order is used instead.

The resulting `BaMM.representation` is a `(5^(order+1), motif_length)` matrix
of `Float32` log-odds scores. Row indexing follows the 5-ary code convention:
`code = b[1] * 5^order + b[2] * 5^(order-1) + ... + b[order+1] * 5^0`.
"""
function read_bamm(path::AbstractString; order::Union{Integer,Nothing}=nothing)
    isfile(path) || throw(ModelFormatError(path, "file not found."))

    raw_blocks = _parse_bamm_blocks(path)
    n_positions = length(raw_blocks)
    n_positions > 0 || throw(ModelFormatError(path, "no valid position blocks found."))

    max_order = length(raw_blocks[1]) - 1
    max_order >= 0 || throw(ModelFormatError(path, "file has no order data."))

    target_order = order === nothing ? max_order : min(Int(order), max_order)
    target_order < 0 &&
        throw(ModelFormatError(path, "order must be non-negative, got $order."))

    # Build the 5-ary tensor representation for each position, then flatten.
    rep = _build_bamm_representation(raw_blocks, target_order, n_positions, path)

    name = _basename_without_extension(path)

    return BaMM(name, rep, target_order, n_positions)
end

"""
    _parse_bamm_blocks(path)

Parse a BaMM file into a vector of position blocks.
Each block is a vector of arrays (one array per order k, containing 4^(k+1) values).
"""
function _parse_bamm_blocks(path::AbstractString)
    blocks = Vector{Vector{Vector{Float32}}}()
    current_block = Vector{Vector{Float32}}()

    open(path, "r") do io
        while !eof(io)
            line = readline(io)
            if length(line) > MAX_LINE_LENGTH
                throw(ModelFormatError(path, "line exceeds length limit."))
            end
            stripped = strip(line)
            if isempty(stripped)
                # Block separator
                if !isempty(current_block)
                    push!(blocks, current_block)
                    current_block = Vector{Vector{Float32}}()
                end
                continue
            end
            if startswith(stripped, '#')
                continue
            end
            # Parse the line as floats
            parts = split(stripped)
            isempty(parts) && continue
            row = Vector{Float32}(undef, length(parts))
            for (j, p) in enumerate(parts)
                v = tryparse(Float32, p)
                v === nothing && throw(ModelFormatError(path, "non-numeric value: $p."))
                row[j] = v
            end
            if !all(isfinite, row)
                throw(ModelFormatError(path, "non-finite values in BaMM data."))
            end
            push!(current_block, row)
        end
    end

    # Don't forget the last block if file doesn't end with blank line
    if !isempty(current_block)
        push!(blocks, current_block)
    end

    # Validate consistency: all blocks must have the same number of order lines
    if !isempty(blocks)
        n_orders = length(blocks[1])
        for (i, block) in enumerate(blocks)
            if length(block) != n_orders
                throw(
                    ModelFormatError(
                        path,
                        "inconsistent orders in block $i: expected $n_orders, got $(length(block)).",
                    ),
                )
            end
        end
    end

    # Validate order widths
    for (pos_idx, block) in enumerate(blocks)
        for (k, arr) in enumerate(block)
            expected_width = 4^k
            if length(arr) != expected_width
                throw(
                    ModelFormatError(
                        path,
                        "BaMM order $(k - 1) width in block $(pos_idx): expected $expected_width, got $(length(arr)).",
                    ),
                )
            end
        end
    end

    n_positions = length(blocks)
    n_positions > MAX_BAMM_POSITIONS && throw(
        ModelFormatError(
            path, "BaMM has $n_positions positions, exceeds limit $MAX_BAMM_POSITIONS."
        ),
    )

    return blocks
end

"""
    _build_bamm_representation(blocks, target_order, n_positions, path)

Build the flattened 2D log-odds representation from parsed blocks.

For each position `pos` (0-indexed):
1. Determine effective order `current_k = min(pos, target_order)`
2. Get probabilities for order `current_k`
3. Compute log-odds: `log((p + eps) / (0.25^(current_k+1) + eps))`
4. Reshape to 4^(current_k+1) and broadcast to 5^(target_order+1) if needed
5. Set N-state entries (any dimension = 4) to per-position minimum
6. Flatten to a row vector of length 5^(target_order+1)
"""
function _build_bamm_representation(
    blocks::Vector{Vector{Vector{Float32}}},
    target_order::Int,
    n_positions::Int,
    path::AbstractString,
)
    n_rows = 5^(target_order + 1)
    rep = Matrix{Float32}(undef, n_rows, n_positions)

    for pos in 0:(n_positions - 1)
        current_k = min(pos, target_order)

        # Get the probability array for this position at order current_k
        # blocks is 1-indexed: blocks[pos+1][current_k+1]
        p_motif = blocks[pos + 1][current_k + 1]
        expected_width = 4^(current_k + 1)
        if length(p_motif) != expected_width
            throw(
                ModelFormatError(
                    path,
                    "position $pos order $current_k: expected $expected_width values, got $(length(p_motif)).",
                ),
            )
        end

        # Compute log-odds against uniform background
        uniform_bg = Float32(0.25)^(current_k + 1)
        log_odds = Vector{Float32}(undef, expected_width)
        @inbounds for i in 1:expected_width
            log_odds[i] = log((p_motif[i] + BAMM_EPSILON) / (uniform_bg + BAMM_EPSILON))
        end

        # Fill the column: first set all entries, then handle N-state
        # ACGT entries: codes where all 5-ary digits are 0-3
        # N entries: codes containing at least one digit = 4
        missing_dims = target_order - current_k
        col_min = minimum(log_odds)  # min over ACGT scores

        @inbounds for row in 1:n_rows
            code = row - 1  # 0-indexed code
            digits = _decode_5ary(code, target_order + 1)
            if any(d -> d == 4, digits)
                # N-state: use per-position minimum
                rep[row, pos + 1] = col_min
            else
                # ACGT entry: use broadcasted log-odds
                # The leading missing_dims digits are broadcast (ignored)
                relevant_digits = digits[(missing_dims + 1):end]
                # Convert to 4-ary index
                idx = 0
                for d in relevant_digits
                    idx = idx * 4 + d
                end
                rep[row, pos + 1] = log_odds[idx + 1]
            end
        end
    end

    return rep
end

"""
    _decode_5ary(code, n_digits)

Decode an integer code into its 5-ary digits (MSB first).
Returns a vector of `n_digits` integers, each in 0..4.
"""
function _decode_5ary(code::Int, n_digits::Int)
    digits = Vector{Int}(undef, n_digits)
    remaining = code
    for i in n_digits:-1:1
        digits[i] = remaining % 5
        remaining = div(remaining, 5)
    end
    return digits
end
