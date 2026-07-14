# Slim XML reader: parses Jstacs Slim (GenDisMix) model files and
# materializes the mixture parameters (component/ancestor/dependency) into a
# dense 5-ary tensor, then flattens to a 2D matrix suitable for the scanning
# kernel.
#
# XML path to model:
#   .//SLIM
#
# Inside SLIM:
#   <length>           — motif length
#   <distance>         — distance parameter (parsed, not used for scoring)
#   <componentMixtureParameters>  — per-position component weights
#   <ancestorMixtureParameters>  — per-position/component ancestor weights
#   <dependencyParameters>       — per-position/component context-conditioned
#                                  log-probability tables
#
# The dense tensor is built exactly as in Python `read_slim`:
#   1. Normalize each component/ancestor/dependency vector with log-sum-exp.
#   2. Compute the span = max(component_index + ancestor_count - 1) over all
#      positions and components with component_index >= 1.
#   3. For each motif position and each concrete ACGT context (4^span),
#      evaluate the per-symbol log-probability via log-sum-exp over
#      components and ancestors (matching `_slim_symbol_log_probs`).
#   4. Materialize the dense `(5,)*（span+1)` column via the shared
#      `_build_position_column` (which also fills N-states with per-axis
#      minima), identical to the Dimont path.
#
# Shared helpers (`_iter_contexts`, `_context_value`, `_build_position_column`,
# `_decode_5ary`, `_encode_5ary`, `_xml_numeric`, `LOG_UNIFORM_BASE`,
# `_basename_without_extension`) are defined in `dimont_reader.jl` /
# `bamm_reader.jl` and reused here.

const SLIM_MAX_LENGTH = 10_000
const SLIM_MAX_SPAN = 10

"""
    read_slim(path)

Read a Jstacs Slim XML model file and return a [`Slim`](@ref).

The XML mixture parameters are materialized into a dense `(5^(span+1),
motif_length)` matrix of `Float32` log-odds scores.
"""
function read_slim(path::AbstractString)
    isfile(path) || throw(ModelFormatError(path, "file not found."))

    root = xml_parse(path)

    slim = xml_find(root, "//SLIM")
    slim === nothing &&
        throw(ModelFormatError(path, "could not find Slim motif model (SLIM)."))

    length_val = Int(_xml_numeric(xml_find(slim, "length"), path, "Slim length"))
    length_val > SLIM_MAX_LENGTH && throw(
        ModelFormatError(path, "Slim length $length_val exceeds limit $SLIM_MAX_LENGTH."),
    )

    # Parse the three parameter arrays.
    component_raw = _slim_parse_component_params(slim, path)
    ancestor_raw = _slim_parse_ancestor_params(slim, path)
    dependency_raw = _slim_parse_dependency_params(slim, path)

    if length(component_raw) != length_val
        throw(
            ModelFormatError(
                path,
                "Slim componentMixtureParameters length ($(length(component_raw))) does not match length ($length_val).",
            ),
        )
    end

    # Compute span.
    span = _slim_span(component_raw, ancestor_raw, path)
    span > SLIM_MAX_SPAN &&
        throw(ModelFormatError(path, "Slim span $span exceeds limit $SLIM_MAX_SPAN."))

    # Normalize parameters to log-probabilities.
    params = _SlimParams(
        component_raw, ancestor_raw, dependency_raw, span, length_val, path
    )

    rep = _build_slim_representation(params, path)
    name = _basename_without_extension(path)
    return Slim(name, rep, span, length_val)
end

# ── Parameter parsing ───────────────────────────────────────────────────

# Return the `<pos>` children of `elem` in document order.
function _slim_pos_children(elem::XMLElement)
    return [c for c in elem.children if c.tag == "pos"]
end

# componentMixtureParameters: pos[position] -> pos[component] -> numeric
# Returns Vector{Vector{Float64}} indexed [position][component].
function _slim_parse_component_params(slim::XMLElement, path::AbstractString)
    elem = xml_find(slim, "componentMixtureParameters")
    elem === nothing && throw(
        ModelFormatError(path, "malformed Slim model: missing componentMixtureParameters."),
    )
    result = Vector{Vector{Float64}}()
    for pos_p in _slim_pos_children(elem)
        comps = _slim_pos_children(pos_p)
        isempty(comps) &&
            throw(ModelFormatError(path, "malformed Slim model: empty component mixture."))
        vals = Float64[_xml_numeric(c, path, "Slim component weight") for c in comps]
        push!(result, vals)
    end
    isempty(result) &&
        throw(ModelFormatError(path, "malformed Slim model: no component positions."))
    return result
end

# ancestorMixtureParameters: pos[position] -> pos[component] -> pos[ancestor] -> numeric
# Returns Vector{Vector{Vector{Float64}}} indexed [position][component][ancestor].
function _slim_parse_ancestor_params(slim::XMLElement, path::AbstractString)
    elem = xml_find(slim, "ancestorMixtureParameters")
    elem === nothing && throw(
        ModelFormatError(path, "malformed Slim model: missing ancestorMixtureParameters."),
    )
    result = Vector{Vector{Vector{Float64}}}()
    for pos_p in _slim_pos_children(elem)
        comps = _slim_pos_children(pos_p)
        comp_list = Vector{Vector{Float64}}()
        for comp_c in comps
            ancs = _slim_pos_children(comp_c)
            isempty(ancs) && throw(
                ModelFormatError(path, "malformed Slim model: empty ancestor mixture.")
            )
            vals = Float64[_xml_numeric(a, path, "Slim ancestor weight") for a in ancs]
            push!(comp_list, vals)
        end
        push!(result, comp_list)
    end
    return result
end

# dependencyParameters: pos[position] -> pos[component] -> pos[row] -> pos[symbol] -> numeric
# Returns Vector{Vector{Vector{Vector{Float64}}}} indexed [position][component][row][symbol].
function _slim_parse_dependency_params(slim::XMLElement, path::AbstractString)
    elem = xml_find(slim, "dependencyParameters")
    elem === nothing &&
        throw(ModelFormatError(path, "malformed Slim model: missing dependencyParameters."))
    result = Vector{Vector{Vector{Vector{Float64}}}}()
    for pos_p in _slim_pos_children(elem)
        comps = _slim_pos_children(pos_p)
        comp_list = Vector{Vector{Vector{Float64}}}()
        for comp_c in comps
            rows = _slim_pos_children(comp_c)
            isempty(rows) && throw(
                ModelFormatError(path, "malformed Slim model: empty dependency rows.")
            )
            row_list = Vector{Vector{Float64}}()
            for row_r in rows
                syms = _slim_pos_children(row_r)
                isempty(syms) && throw(
                    ModelFormatError(path, "malformed Slim model: empty dependency row."),
                )
                vals = Float64[_xml_numeric(s, path, "Slim dependency value") for s in syms]
                push!(row_list, vals)
            end
            push!(comp_list, row_list)
        end
        push!(result, comp_list)
    end
    return result
end

# ── Span computation ────────────────────────────────────────────────────

# span = max over positions, components (component_index >= 1) of
#   (component_index + ancestor_count - 1)
# where component_index is 0-based and ancestor_count = number of ancestors
# for that component. Mirrors Python `_slim_span`.
function _slim_span(
    component_raw::Vector{Vector{Float64}},
    ancestor_raw::Vector{Vector{Vector{Float64}}},
    path::AbstractString,
)
    span = 0
    for position in 1:length(component_raw)
        n_components = length(component_raw[position])
        # component_index is 0-based; Python iterates range(1, len).
        for component_index in 1:(n_components - 1)
            if component_index > length(ancestor_raw[position])
                throw(
                    ModelFormatError(
                        path,
                        "malformed Slim model: missing ancestor mixture for position $position component $component_index.",
                    ),
                )
            end
            ancestor_count = length(ancestor_raw[position][component_index + 1])
            if ancestor_count <= 0
                throw(
                    ModelFormatError(
                        path,
                        "malformed Slim model: empty ancestor mixture at position $(position - 1).",
                    ),
                )
            end
            # 0-based component_index + ancestor_count - 1
            reach = (component_index) + ancestor_count - 1
            span = max(span, reach)
        end
    end
    return span
end

# ── Normalized parameters ───────────────────────────────────────────────

struct _SlimParams
    # component_log_probs[position] :: Vector{Float64}
    component_log_probs::Vector{Vector{Float64}}
    # ancestor_log_probs[position][component] :: Vector{Float64}
    ancestor_log_probs::Vector{Vector{Vector{Float64}}}
    # dependency_log_probs[position][component] :: Matrix{Float64} (rows=contexts, cols=symbols)
    dependency_log_probs::Vector{Vector{Matrix{Float64}}}
    span::Int
    length::Int
    # alphabet_size (= number of symbols per dependency row, always 4 for DNA)
    alphabet_size::Int
end

function _SlimParams(
    component_raw::Vector{Vector{Float64}},
    ancestor_raw::Vector{Vector{Vector{Float64}}},
    dependency_raw::Vector{Vector{Vector{Vector{Float64}}}},
    span::Int,
    length_val::Int,
    path::AbstractString,
)
    # alphabet_size = number of symbols in the first dependency row.
    if isempty(dependency_raw) ||
        isempty(dependency_raw[1]) ||
        isempty(dependency_raw[1][1]) ||
        isempty(dependency_raw[1][1][1])
        throw(ModelFormatError(path, "malformed Slim model: empty dependency parameters."))
    end
    alphabet_size = length(dependency_raw[1][1][1])

    component_log_probs = [_log_normalize(component_raw[p]) for p in 1:length_val]
    ancestor_log_probs = [
        [_log_normalize(ancestor_raw[p][c]) for c in 1:length(ancestor_raw[p])] for
        p in 1:length_val
    ]
    dependency_log_probs = [
        [_log_normalize_rows(dependency_raw[p][c]) for c in 1:length(dependency_raw[p])] for
        p in 1:length_val
    ]

    return _SlimParams(
        component_log_probs,
        ancestor_log_probs,
        dependency_log_probs,
        span,
        length_val,
        alphabet_size,
    )
end

# Normalize each row of a dependency table (rows = context codes, cols = symbols).
function _log_normalize_rows(rows::Vector{Vector{Float64}})
    n_rows = length(rows)
    n_cols = length(rows[1])
    mat = Matrix{Float64}(undef, n_rows, n_cols)
    for r in 1:n_rows
        if length(rows[r]) != n_cols
            throw(ModelFormatError("", "inconsistent dependency row lengths."))
        end
        normalized = _log_normalize(rows[r])
        for c in 1:n_cols
            mat[r, c] = normalized[c]
        end
    end
    return mat
end

# ── log-sum-exp helpers (Float64) ───────────────────────────────────────

function _logsumexp(v::AbstractVector{<:Real})
    m = maximum(v)
    if !isfinite(m)
        # all -Inf: return -Inf; if any +Inf, m=+Inf and result is +Inf
        return m
    end
    s = 0.0
    for x in v
        s += exp(x - m)
    end
    return m + log(s)
end

function _log_normalize(v::AbstractVector{<:Real})
    lse = _logsumexp(v)
    return Float64[x - lse for x in v]
end

# ── Per-symbol log-probability (mirrors `_slim_symbol_log_probs`) ────────

# position is 0-based, matching Python. full_context is a tuple of length span
# of concrete ACGT indices (0..3).
function _slim_symbol_log_probs(
    position::Int,
    symbol::Int,
    full_context::Tuple,
    params::_SlimParams,
    path::AbstractString,
)
    comp_lp = params.component_log_probs[position + 1]
    dep_lp = params.dependency_log_probs[position + 1]
    anc_lp = params.ancestor_log_probs[position + 1]
    span = params.span
    alphabet = params.alphabet_size

    n_comp = length(comp_lp)
    local_scores = Vector{Float64}(undef, n_comp)

    # Component 0 (0-based): no ancestors, context row 0.
    local_scores[1] = comp_lp[1] + dep_lp[1][1, symbol + 1]

    for component_index in 1:(n_comp - 1)
        ancestor_count = length(anc_lp[component_index + 1])
        context_index = 0
        # current_order in 1..(component_index-1) (0-based component_index)
        for current_order in 1:(component_index - 1)
            parent_position = position - current_order
            ctx_val = _context_value(full_context, span, position, parent_position, path)
            context_index = context_index * alphabet + ctx_val
        end

        dependency = dep_lp[component_index + 1]
        total_contexts = size(dependency, 1)
        width = size(dependency, 2)

        ancestor_scores = Vector{Float64}(undef, ancestor_count)
        for ancestor_index in 0:(ancestor_count - 1)
            parent_position = position - component_index - ancestor_index
            ctx_val = _context_value(full_context, span, position, parent_position, path)
            context_index = mod1(context_index * width + ctx_val + 1, total_contexts) - 1
            ancestor_scores[ancestor_index + 1] =
                anc_lp[component_index + 1][ancestor_index + 1] +
                dependency[context_index + 1, symbol + 1]
        end
        local_scores[component_index + 1] =
            comp_lp[component_index + 1] + _logsumexp(ancestor_scores)
    end

    return _logsumexp(local_scores) + LOG_UNIFORM_BASE
end

# ── Dense representation materialization ─────────────────────────────────

function _build_slim_representation(params::_SlimParams, path::AbstractString)
    span = params.span
    length_val = params.length
    n_rows = 5^(span + 1)
    rep = Matrix{Float32}(undef, n_rows, length_val)

    full_contexts = _iter_contexts(span)

    for position in 0:(length_val - 1)
        # context_scores[ctx] = Vector{Float64} of 4 symbol log-probs
        context_scores = Dict{Tuple{Vararg{Int}},Vector{Float64}}()
        for ctx in full_contexts
            symbol_log_probs = Vector{Float64}(undef, 4)
            for symbol in 0:3
                symbol_log_probs[symbol + 1] = _slim_symbol_log_probs(
                    position, symbol, ctx, params, path
                )
            end
            context_scores[ctx] = symbol_log_probs
        end
        _build_position_column(context_scores, span, position + 1, rep)
    end

    return rep
end
