# Score profile reader: FASTA-like numerical score files.

"""
    read_scores(path::AbstractString)

Read a FASTA-like file of numerical score profiles into a [`ScoreProfile`](@ref).

Each `>` header starts a new profile. Subsequent lines contain
whitespace- or comma-separated float values. The model name is derived
from the filename (without extension).
"""
const MAX_SCORE_FILE_BYTES = 256 * 1024^2
const MAX_SCORE_LINE_BYTES = 4 * 1024^2
const MAX_SCORE_ROWS = 1_000_000
const MAX_SCORE_ELEMENTS = 100_000_000

function read_scores(path::AbstractString)
    file = String(path)
    isfile(file) || throw(ModelFormatError(file, "score file does not exist."))
    filesize(file) <= MAX_SCORE_FILE_BYTES ||
        throw(ModelFormatError(file, "score file exceeds the size limit."))

    rows = Vector{Vector{Float32}}()
    current_values = Float32[]
    seen_header = false
    elements = 0
    try
        for line in eachline(file)
            ncodeunits(line) <= MAX_SCORE_LINE_BYTES ||
                throw(ModelFormatError(file, "score line exceeds the size limit."))
            stripped = strip(line)
            isempty(stripped) && continue
            if startswith(stripped, '>')
                seen_header && push!(rows, current_values)
                length(rows) < MAX_SCORE_ROWS ||
                    throw(ModelFormatError(file, "score row count exceeds the limit."))
                current_values = Float32[]
                seen_header = true
                continue
            end
            seen_header || throw(ModelFormatError(file, "score values require a header."))
            cleaned = replace(stripped, ',' => ' ')
            for token in split(cleaned)
                isempty(token) && continue
                value = try
                    parse(Float32, String(token))
                catch
                    throw(ModelFormatError(file, "invalid score value: '$token'."))
                end
                isfinite(value) ||
                    throw(ModelFormatError(file, "score values must be finite."))
                elements += 1
                elements <= MAX_SCORE_ELEMENTS ||
                    throw(ModelFormatError(file, "score element count exceeds the limit."))
                push!(current_values, value)
            end
        end
    catch error
        error isa ModelFormatError && rethrow()
        throw(
            ModelFormatError(
                file, "could not read score file: $(sprint(showerror, error))."
            ),
        )
    end
    seen_header && push!(rows, current_values)
    length(rows) <= MAX_SCORE_ROWS ||
        throw(ModelFormatError(file, "score row count exceeds the limit."))
    isempty(rows) && throw(ModelFormatError(file, "score file contains no profiles."))

    name = splitext(basename(file))[1]
    return ScoreProfile(name, build_ragged(rows))
end

@doc "Read a bounded FASTA-like numerical score profile file." read_scores
