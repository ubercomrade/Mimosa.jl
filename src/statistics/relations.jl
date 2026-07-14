# Group relation parsing for null-distribution builds.
#
# Reads a TSV/CSV file with motif name and group columns, builds a mapping
# from each motif to the set of motifs from a *different* group.

"""
    GroupRelations

Mapping from motif name to the set of eligible comparison targets (motifs
from a different group).

Fields:
- `groups::Dict{String,String}`: motif name → group name.
- `eligible::Dict{String,Set{String}}`: motif name → set of eligible targets.
"""
struct GroupRelations
    groups::Dict{String,String}
    eligible::Dict{String,Set{String}}
end

"""
    parse_group_relations(path; name_column="motif", group_column="group", ignore_missing=false, known_names=nothing)

Parse a motif-to-group table (TSV/CSV) and build eligible pairs: for each
motif, all other motifs from a *different* group.

# Arguments
- `path`: path to the relation file.
- `name_column`: column name for motif names (default `"motif"`).
- `group_column`: column name for group names (default `"group"`).
- `ignore_missing`: if `true`, silently skip motif names not in `known_names`.
- `known_names`: optional set of known motif names for validation.

Returns a [`GroupRelations`](@ref) with the group mapping and eligible pairs.
"""
function parse_group_relations(
    path::AbstractString;
    name_column::AbstractString="motif",
    group_column::AbstractString="group",
    ignore_missing::Bool=false,
    known_names::Union{Nothing,Set{String}}=nothing,
)
    rows, headers = _read_table(path)

    name_idx = findfirst(==(name_column), headers)
    group_idx = findfirst(==(group_column), headers)

    if name_idx === nothing
        throw(
            ArgumentError(
                "Group table must contain column '$name_column'. Found: $(join(headers, ", "))",
            ),
        )
    end
    if group_idx === nothing
        throw(
            ArgumentError(
                "Group table must contain column '$group_column'. Found: $(join(headers, ", "))",
            ),
        )
    end

    groups = Dict{String,String}()
    for row in rows
        name = String(strip(row[name_idx]))
        group = String(strip(row[group_idx]))
        isempty(name) && continue
        isempty(group) && throw(ArgumentError("Relation group names must not be empty."))
        haskey(groups, name) &&
            throw(ArgumentError("Relation file contains duplicate motif '$name'."))
        groups[name] = group
    end

    # Validate against known names
    if known_names !== nothing
        missing_names = setdiff(keys(groups), known_names)
        if !isempty(missing_names) && !ignore_missing
            throw(
                ArgumentError(
                    "Relation input references unknown motifs: $(join(sort!(collect(missing_names)), ", "))",
                ),
            )
        end
    end

    # Build eligible pairs
    names = sort!(collect(keys(groups)))
    if known_names !== nothing
        names = filter(n -> n in known_names, names)
    end

    eligible = Dict{String,Set{String}}()
    for query in names
        targets = Set{String}()
        for target in names
            if target != query && groups[target] != groups[query]
                push!(targets, target)
            end
        end
        eligible[query] = targets
    end

    return GroupRelations(groups, eligible)
end

"""
    eligible_targets(relations::GroupRelations, query::AbstractString)

Return the sorted list of eligible target names for a given query motif.
"""
function eligible_targets(relations::GroupRelations, query::AbstractString)
    targets = get(relations.eligible, query, Set{String}())
    return sort!(collect(targets))
end

# ---------------------------------------------------------------------------
# Table reader (TSV/CSV with delimiter sniffing)
# ---------------------------------------------------------------------------

function _read_table(path::AbstractString)
    isfile(path) ||
        throw(ArgumentError("could not read relation file '$path': file not found."))
    filesize(path) <= 256 * 1024^2 ||
        throw(ArgumentError("relation file exceeds the size limit."))
    sample = open(path, "r") do io
        return read(io, min(filesize(path), 4096))
    end
    delimiter = _sniff_delimiter(String(sample))
    lines = String[]
    open(path, "r") do io
        for line in eachline(io)
            ncodeunits(line) <= 4 * 1024^2 ||
                throw(ArgumentError("relation line exceeds the size limit."))
            stripped = strip(line)
            isempty(stripped) || startswith(stripped, '#') || push!(lines, line)
        end
    end
    isempty(lines) && throw(ArgumentError("Relation file is empty: $path"))

    headers = String.(strip.(split(lines[1], delimiter; keepempty=true)))
    isempty(headers) && throw(ArgumentError("Relation file has no header: $path"))
    any(isempty, headers) &&
        throw(ArgumentError("Relation file contains an empty header: $path"))
    length(unique(headers)) == length(headers) ||
        throw(ArgumentError("Relation file contains duplicate headers: $path"))

    rows = Vector{Vector{String}}()
    for line in lines[2:end]
        fields = split(line, delimiter; keepempty=true)
        length(fields) == length(headers) || throw(
            ArgumentError(
                "Relation row has $(length(fields)) fields; expected $(length(headers)).",
            ),
        )
        row = String.(strip.(fields))
        push!(rows, row)
    end

    return rows, headers
end

function _sniff_delimiter(content::AbstractString)
    # Try to detect delimiter: tab, comma, or semicolon
    sample = first(content, min(length(content), 4096))
    tab_count = count(==('\t'), sample)
    comma_count = count(==(','), sample)
    semicolon_count = count(==(';'), sample)

    if tab_count >= comma_count && tab_count >= semicolon_count && tab_count > 0
        return '\t'
    elseif comma_count > 0 && comma_count >= semicolon_count
        return ','
    elseif semicolon_count > 0
        return ';'
    else
        return '\t'
    end
end
