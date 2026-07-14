# FASTA reader producing EncodedSequenceBatch.

const MAX_FASTA_SEQUENCES = 1_000_000
const MAX_FASTA_SEQUENCE_LENGTH = 100_000_000
const MAX_FASTA_LINE_LENGTH = 1_000_000

"""
    read_fasta(path; max_sequences=MAX_FASTA_SEQUENCES)

Read a FASTA file and return `(batch, names)` where `batch` is an
[`EncodedSequenceBatch`](@ref) and `names` is a `Vector{String}` of sequence
identifiers (first whitespace-delimited token after `>`).

Sequences are encoded using 5-ary encoding (A=0, C=1, G=2, T=3, N/ambiguous=4).
Lowercase is normalized to uppercase before encoding. All non-ACGT characters
map to 4, matching Python behavior exactly.

Empty sequences are allowed.
"""
function read_fasta(path::AbstractString; max_sequences::Int=MAX_FASTA_SEQUENCES)
    isfile(path) || throw(ModelFormatError(path, "file not found."))
    return open(path, "r") do io
        return _read_fasta_io(io, path, max_sequences)
    end
end

"""
    readsequences(path; kwargs...)

Alias for [`read_fasta`](@ref). Reads a FASTA file and returns
`(batch, names)`.
"""
readsequences(path::AbstractString; kwargs...) = read_fasta(path; kwargs...)

function _read_fasta_io(io::IO, path::AbstractString, max_sequences::Int)
    rows = Vector{Vector{UInt8}}()
    names = Vector{String}()
    current_seq = UInt8[]
    current_name = ""
    has_current = false
    n_sequences = 0

    while !eof(io)
        line = readline(io)
        if length(line) > MAX_FASTA_LINE_LENGTH
            throw(
                ModelFormatError(path, "line exceeds length limit $MAX_FASTA_LINE_LENGTH.")
            )
        end
        stripped = strip(line)
        isempty(stripped) && continue

        if stripped[1] == '>'
            if has_current
                push!(rows, current_seq)
                push!(names, current_name)
                n_sequences += 1
                n_sequences >= max_sequences && throw(
                    ModelFormatError(path, "exceeded max_sequences limit $max_sequences."),
                )
                current_seq = UInt8[]
            end
            # Header: take the first whitespace-delimited token as the name
            header = strip(stripped[2:end])
            if isempty(header)
                current_name = ""
            else
                name_parts = split(header, isspace; limit=2)
                current_name = isempty(name_parts) ? "" : String(name_parts[1])
            end
            has_current = true
        else
            has_current ||
                throw(ModelFormatError(path, "sequence data before header line."))
            _append_encoded!(current_seq, stripped)
            length(current_seq) > MAX_FASTA_SEQUENCE_LENGTH && throw(
                ModelFormatError(
                    path, "sequence exceeds length limit $MAX_FASTA_SEQUENCE_LENGTH."
                ),
            )
        end
    end

    if has_current
        push!(rows, current_seq)
        push!(names, current_name)
    end

    isempty(rows) && throw(ModelFormatError(path, "no sequences found in FASTA file."))

    return (EncodedSequenceBatch(rows), names)
end

# Append ASCII-encoded nucleotides from a line to the current sequence buffer.
function _append_encoded!(seq::Vector{UInt8}, line::AbstractString)
    for i in 1:ncodeunits(line)
        push!(seq, _ENCODE_TABLE[codeunit(line, i) + 1])
    end
    return seq
end
