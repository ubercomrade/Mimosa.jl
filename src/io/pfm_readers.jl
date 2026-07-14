# MEME and PFM parsers with size limits and clear errors.

const MAX_MEME_MOTIF_LENGTH = 10_000
const MAX_PFM_LENGTH = 10_000
const MAX_LINE_LENGTH = 1_000_000

"""
    read_meme(path; index=0)

Read one motif from a MEME letter-probability file and return a frequency matrix.

`index` selects the motif by file order (zero-based for compatibility with
the Python API; internally converted to one-based).
"""
function read_meme(path::AbstractString; index::Integer=0)
    isfile(path) || throw(ModelFormatError(path, "file not found."))
    idx = Int(index)
    idx < 0 &&
        throw(ModelFormatError(path, "motif index must be non-negative, got $index."))
    return open(path, "r") do io
        return _read_meme_io(io, path, idx)
    end
end

function _read_meme_io(io::IO, path::AbstractString, target_index::Int)
    motif_count = 0
    while !eof(io)
        line = readline(io)
        if startswith(line, "MOTIF")
            is_target = motif_count == target_index
            motif_count += 1
            parts = split(strip(line))
            if length(parts) < 2
                throw(ModelFormatError(path, "MOTIF line has no name."))
            end
            name = parts[2]
            motif_length = _meme_length_from_header(readline(io), path, name)
            if is_target
                matrix = _read_meme_matrix_rows(io, path, name, motif_length)
                _validate_meme_matrix(matrix, path, name, motif_length)
                pfm = _orient_meme_matrix(matrix, motif_length)
                return (name=String(name), frequencies=pfm)
            else
                motif_length <= 0 &&
                    throw(ModelFormatError(path, "motif $(name) has invalid length."))
                for _ in 1:motif_length
                    eof(io) && throw(
                        ModelFormatError(
                            path, "motif $(name) has fewer rows than declared length."
                        ),
                    )
                    readline(io)
                end
            end
        end
    end
    if motif_count == 0
        throw(ModelFormatError(path, "no motifs found."))
    end
    return throw(
        ModelFormatError(
            path,
            "motif index $(target_index) out of range. File contains $(motif_count) motifs.",
        ),
    )
end

function _meme_length_from_header(
    header_line::AbstractString, path::AbstractString, name::AbstractString
)
    header = split(strip(header_line))
    idx = findfirst(==("w="), header)
    if idx === nothing || idx >= length(header)
        throw(
            ModelFormatError(path, "motif $(name) header has no valid 'w=' length field.")
        )
    end
    motif_length = tryparse(Int, header[idx + 1])
    if motif_length === nothing
        throw(
            ModelFormatError(
                path, "motif $(name) length is not an integer: $(header[idx + 1])."
            ),
        )
    end
    return motif_length
end

function _read_meme_matrix_rows(
    io::IO, path::AbstractString, name::AbstractString, nrows::Int
)
    nrows <= 0 && throw(ModelFormatError(path, "motif $(name) has invalid length."))
    nrows > MAX_MEME_MOTIF_LENGTH && throw(
        ModelFormatError(
            path, "motif $(name) length $nrows exceeds limit $MAX_MEME_MOTIF_LENGTH."
        ),
    )
    rows = Vector{Vector{Float32}}(undef, nrows)
    for i in 1:nrows
        eof(io) && throw(
            ModelFormatError(path, "motif $(name) has fewer rows than declared length.")
        )
        line = readline(io)
        if length(line) > MAX_LINE_LENGTH
            throw(ModelFormatError(path, "motif $(name) row exceeds line length limit."))
        end
        parts = split(strip(line))
        if length(parts) != NUCLEOTIDE_CARDINALITY
            throw(
                ModelFormatError(
                    path,
                    "motif $(name) row $i has $(length(parts)) columns, expected $NUCLEOTIDE_CARDINALITY.",
                ),
            )
        end
        row = Vector{Float32}(undef, NUCLEOTIDE_CARDINALITY)
        for (j, p) in enumerate(parts)
            v = tryparse(Float32, p)
            v === nothing && throw(
                ModelFormatError(path, "motif $(name) row $i has non-numeric value: $(p)."),
            )
            row[j] = v
        end
        rows[i] = row
    end
    return rows
end

function _validate_meme_matrix(
    rows::Vector{Vector{Float32}}, path::AbstractString, name::AbstractString, nrows::Int
)
    for (i, row) in enumerate(rows)
        if !all(isfinite, row)
            throw(
                ModelFormatError(
                    path, "motif $(name) contains non-finite values in row $i."
                ),
            )
        end
    end
    return nothing
end

# MEME stores `position × base`; transpose to `(base, position)` PFM.
function _orient_meme_matrix(rows::Vector{Vector{Float32}}, npos::Int)
    pfm = Matrix{Float32}(undef, NUCLEOTIDE_CARDINALITY, npos)
    for (pos, row) in enumerate(rows)
        for base in 1:NUCLEOTIDE_CARDINALITY
            pfm[base, pos] = row[base]
        end
    end
    return pfm
end

"""
    read_pfm(path)

Read a Position Frequency Matrix from a plain-text file.

Auto-detects orientation: if the file has 4 or 5 columns it is transposed to
`(base, position)`; if it has 4 or 5 rows it is kept as-is.
"""
function read_pfm(path::AbstractString)
    isfile(path) || throw(ModelFormatError(path, "file not found."))
    filesize(path) <= 256 * 1024^2 ||
        throw(ModelFormatError(path, "PFM file exceeds the size limit."))
    rows = Vector{Vector{Float32}}()
    open(path, "r") do io
        while !eof(io)
            line = readline(io)
            length(line) > MAX_LINE_LENGTH &&
                throw(ModelFormatError(path, "line exceeds length limit."))
            stripped = strip(line)
            isempty(stripped) && continue
            startswith(stripped, ">") && continue
            parts = split(stripped)
            row = Vector{Float32}(undef, length(parts))
            for (j, p) in enumerate(parts)
                v = tryparse(Float32, p)
                v === nothing && throw(ModelFormatError(path, "non-numeric value: $(p)."))
                row[j] = v
            end
            length(rows) < MAX_PFM_LENGTH ||
                throw(ModelFormatError(path, "PFM row count exceeds the size limit."))
            length(parts) <= MAX_PFM_LENGTH ||
                throw(ModelFormatError(path, "PFM column count exceeds the size limit."))
            push!(rows, row)
        end
    end
    isempty(rows) && throw(ModelFormatError(path, "PFM file is empty."))
    ncols = length(rows[1])
    all(r -> length(r) == ncols, rows) ||
        throw(ModelFormatError(path, "PFM rows have inconsistent column counts."))
    n_rows = length(rows)
    raw = Matrix{Float32}(undef, n_rows, ncols)
    for (i, row) in enumerate(rows)
        for j in 1:ncols
            raw[i, j] = row[j]
        end
    end
    n_rows > MAX_PFM_LENGTH && throw(
        ModelFormatError(path, "PFM dimension $n_rows exceeds limit $MAX_PFM_LENGTH.")
    )
    ncols > MAX_PFM_LENGTH &&
        throw(ModelFormatError(path, "PFM dimension $ncols exceeds limit $MAX_PFM_LENGTH."))
    if ncols in (NUCLEOTIDE_CARDINALITY, 5)
        pfm = Matrix(transpose(raw))
    elseif n_rows in (NUCLEOTIDE_CARDINALITY, 5)
        pfm = copy(raw)
    else
        throw(ModelFormatError(path, "one axis must contain 4 or 5 nucleotide rows."))
    end
    if size(pfm, 1) == 5
        pfm = pfm[1:4, :]
    end
    if size(pfm, 2) <= 0
        throw(ModelFormatError(path, "motif length must be positive."))
    end
    if !all(isfinite, pfm)
        throw(ModelFormatError(path, "matrix contains non-finite values."))
    end
    name = _basename_without_extension(path)
    return (name=name, frequencies=pfm)
end

function _basename_without_extension(path::AbstractString)
    base = basename(path)
    dot = findlast('.', base)
    return dot === nothing ? base : base[1:(dot - 1)]
end
