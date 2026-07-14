# SiteGA (.mat) parser: reads and writes SiteGA dinucleotide motif files.
#
# File format:
#   Line 1: model name
#   Line 2: LPD count (number of segments)
#   Line 3: model length
#   Line 4: minimum score
#   Line 5: maximum score ("Razmah")
#   Lines 6+: segments with 5 tab-separated fields:
#     start  stop  total_value  dinuc_index  dinucleotide
#
# Each segment distributes `total_value / (stop - start + 1)` across positions
# [start, stop] (0-indexed) for the given dinucleotide pair.
#
# The representation is (5, 5, length) in Python, flattened to (25, length) in
# Julia with row indexing `code = nuc1 * 5 + nuc2`.

using Printf

const MAX_SITEGA_LENGTH = 10_000
const SITEGA_EPSILON::Float64 = 1e-9
const DINUC_MAP = Dict(
    "aa" => (0, 0),
    "ac" => (0, 1),
    "ag" => (0, 2),
    "at" => (0, 3),
    "ca" => (1, 0),
    "cc" => (1, 1),
    "cg" => (1, 2),
    "ct" => (1, 3),
    "ga" => (2, 0),
    "gc" => (2, 1),
    "gg" => (2, 2),
    "gt" => (2, 3),
    "ta" => (3, 0),
    "tc" => (3, 1),
    "tg" => (3, 2),
    "tt" => (3, 3),
)

const DINUC_LIST = [
    "aa",
    "ac",
    "ag",
    "at",
    "ca",
    "cc",
    "cg",
    "ct",
    "ga",
    "gc",
    "gg",
    "gt",
    "ta",
    "tc",
    "tg",
    "tt",
]

"""
    read_sitega(path)

Read a SiteGA motif from a `.mat` file and return a [`SiteGA`](@ref).

The representation is a `(25, motif_length)` `Float32` matrix with row indexing
`code = nuc1 * 5 + nuc2`.
"""
function read_sitega(path::AbstractString)
    isfile(path) || throw(ModelFormatError(path, "file not found."))
    filesize(path) <= 256 * 1024^2 ||
        throw(ModelFormatError(path, "SiteGA file exceeds the size limit."))

    lines = readlines(path)

    # Line 1: model name
    isempty(lines) && throw(ModelFormatError(path, "empty file: missing model name."))
    name = strip(lines[1])
    isempty(name) && throw(ModelFormatError(path, "missing model name."))

    # Line 2: LPD count
    length(lines) < 2 && throw(ModelFormatError(path, "missing LPD count line."))
    lpd_parts = split(strip(lines[2]))
    isempty(lpd_parts) && throw(ModelFormatError(path, "missing LPD count."))
    lpd_count = tryparse(Int, lpd_parts[1])
    lpd_count === nothing &&
        throw(ModelFormatError(path, "invalid LPD count: $(lpd_parts[1])."))

    # Line 3: model length
    length(lines) < 3 && throw(ModelFormatError(path, "missing model length line."))
    len_parts = split(strip(lines[3]))
    isempty(len_parts) && throw(ModelFormatError(path, "missing model length."))
    motif_length = tryparse(Int, len_parts[1])
    motif_length === nothing &&
        throw(ModelFormatError(path, "invalid model length: $(len_parts[1])."))
    motif_length <= 0 &&
        throw(ModelFormatError(path, "model length must be positive, got $motif_length."))
    motif_length > MAX_SITEGA_LENGTH && throw(
        ModelFormatError(
            path, "model length $motif_length exceeds limit $MAX_SITEGA_LENGTH."
        ),
    )

    # Lines 4-5: minimum and maximum (skip)
    # Line 4: minimum score
    # Line 5: maximum score ("Razmah")
    # These are derived from the representation; we skip them.

    # Initialize representation: (25, motif_length) filled with zeros
    rep = Matrix{Float32}(undef, 25, motif_length)
    fill!(rep, 0.0f0)

    # Parse segment lines (starting from line 6, index 6 in 1-based)
    for (line_idx, line) in enumerate(lines[6:end])
        stripped = strip(line)
        isempty(stripped) && continue

        parts = split(stripped)
        if length(parts) != 5
            throw(
                ModelFormatError(
                    path,
                    "line $(line_idx + 6) must contain 5 fields, got $(length(parts)).",
                ),
            )
        end

        start_idx = tryparse(Int, parts[1])
        stop_idx = tryparse(Int, parts[2])
        value = tryparse(Float64, parts[3])
        dinucleotide = lowercase(parts[5])

        start_idx === nothing &&
            throw(ModelFormatError(path, "invalid start index: $(parts[1])."))
        stop_idx === nothing &&
            throw(ModelFormatError(path, "invalid stop index: $(parts[2])."))
        value === nothing && throw(ModelFormatError(path, "invalid value: $(parts[3])."))

        if !haskey(DINUC_MAP, dinucleotide)
            throw(
                ModelFormatError(
                    path, "invalid dinucleotide: $(dinucleotide) on line $(line_idx + 6)."
                ),
            )
        end

        if start_idx < 0 || stop_idx < start_idx || stop_idx >= motif_length
            throw(
                ModelFormatError(
                    path,
                    "range $start_idx-$stop_idx is outside model length $motif_length on line $(line_idx + 6).",
                ),
            )
        end

        nuc1, nuc2 = DINUC_MAP[dinucleotide]
        row_code = nuc1 * 5 + nuc2  # 0-indexed
        n_positions = stop_idx - start_idx + 1
        per_position = Float32(value / n_positions)

        for idx in start_idx:stop_idx
            rep[row_code + 1, idx + 1] += per_position
        end
    end

    if !all(isfinite, rep)
        throw(ModelFormatError(path, "representation contains non-finite values."))
    end

    return SiteGA(name, rep, motif_length)
end

"""
    write_sitega(path, model::SiteGA)

Write a SiteGA motif to a `.mat` file.

The writer groups contiguous positions with identical values into segments,
matching the Python writer format. Only dinucleotides with non-zero entries
(within `SITEGA_EPSILON`) are written.
"""
function write_sitega(path::AbstractString, model::SiteGA)
    rep = model.representation
    motif_length = model.motif_length
    mn, mx = scorebounds(model)

    segments = []

    for nuc1 in 0:3
        for nuc2 in 0:3
            row_code = nuc1 * 5 + nuc2  # 0-indexed
            # Check if this dinucleotide has any non-zero entries
            row_data = rep[row_code + 1, :]
            if all(abs.(row_data) .<= SITEGA_EPSILON)
                continue
            end

            dinucleotide = DINUC_LIST[nuc1 * 4 + nuc2 + 1]
            pos = 1  # 1-indexed

            while pos <= motif_length
                # Skip near-zero positions
                while pos <= motif_length && abs(row_data[pos]) <= SITEGA_EPSILON
                    pos += 1
                end
                if pos > motif_length
                    break
                end

                start_pos = pos
                current_val = row_data[pos]

                # Extend while values are approximately equal
                while pos + 1 <= motif_length &&
                      abs(row_data[pos + 1] - current_val) < SITEGA_EPSILON
                    pos += 1
                end

                push!(
                    segments,
                    (
                        start=start_pos - 1,
                        stop=pos - 1,
                        val=current_val,
                        dinuc=dinucleotide,
                    ),
                )

                pos += 1
            end
        end
    end

    lpd_count = length(segments)

    open(path, "w") do io
        println(io, model.name)
        println(io, "$lpd_count\tLPD count")
        println(io, "$motif_length\tModel length")
        @printf(io, "%.12f\tMinimum\n", mn)
        @printf(io, "%.12f\tRazmah\n", mx)

        # dinuc_index map (lowercase dinucleotide → index, matches Python's itertools.product("acgt", repeat=2))
        dinuc_index = Dict(
            "aa" => 0,
            "ac" => 1,
            "ag" => 2,
            "at" => 3,
            "ca" => 4,
            "cc" => 5,
            "cg" => 6,
            "ct" => 7,
            "ga" => 8,
            "gc" => 9,
            "gg" => 10,
            "gt" => 11,
            "ta" => 12,
            "tc" => 13,
            "tg" => 14,
            "tt" => 15,
        )

        for seg in segments
            range_length = seg.stop - seg.start + 1
            total_value = seg.val * range_length
            idx = dinuc_index[seg.dinuc]
            @printf(
                io,
                "%d\t%d\t%.12f\t%d\t%s\n",
                seg.start,
                seg.stop,
                total_value,
                idx,
                seg.dinuc
            )
        end
    end
    return nothing
end
