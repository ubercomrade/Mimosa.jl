# Dimont XML reader: parses Jstacs Dimont model files and materializes
# the Bayesian network parameter trees into a dense 5-ary tensor, then
# flattens to a 2D matrix suitable for the scanning kernel.
#
# XML path to model:
#   .//ThresholdedStrandChIPper/function/pos/MarkovModelDiffSM
#
# Inside MarkovModelDiffSM:
#   bayesianNetworkSF/trees — one <pos val="N"> per motif position, each
#   containing a <parameterTree> with:
#     - <contextPoss>: list of parent positions for this position
#     - <root>/<treeElement>: recursive tree of scoring parameters
#
# Each treeElement has:
#   - <contextPos>: position used for branching at this node
#   - Either <pars> (leaf: 4 <parameter> entries with <symbol> and <value>)
#   - Or <children> (internal: 4 child treeElements indexed by <pos val="N">)
#
# Materialization:
#   For each motif position, for each concrete ACGT context (4^span combos):
#     Walk the tree by looking up context values at contextPos nodes.
#     At the leaf, get 4 scores (one per symbol).
#     Add log(4.0) to each score (log-uniform base).
#   Fill N-states (any dimension = 4) with the minimum over concrete nucleotides.

include("xml_parser.jl")

const DIMONT_MAX_LENGTH = 10_000
const DIMONT_MAX_SPAN = 10
const LOG_UNIFORM_BASE::Float64 = log(4.0)

# Tree node representation for Dimont parameter trees.
# A leaf node has `scores` (Vector{Float64} of length 4).
# An internal node has `context_pos` and `children` (Vector of 4 TreeNodes).
struct DimontTreeNode
    context_pos::Int
    scores::Union{Nothing,Vector{Float64}}
    children::Union{Nothing,Vector{DimontTreeNode}}
end

"""
    read_dimont(path)

Read a Jstacs Dimont XML model file and return a [`Dimont`](@ref).

The XML parameter trees are materialized into a dense `(5^(span+1),
motif_length)` matrix of `Float32` log-odds scores.
"""
function read_dimont(path::AbstractString)
    isfile(path) || throw(ModelFormatError(path, "file not found."))

    root = xml_parse(path)

    # Find the MarkovModelDiffSM element
    model = xml_find(root, "//MarkovModelDiffSM")
    model === nothing && throw(
        ModelFormatError(path, "could not find Dimont motif model (MarkovModelDiffSM).")
    )

    # Find the trees container
    trees = xml_find(model, "bayesianNetworkSF/trees")
    trees === nothing && throw(
        ModelFormatError(path, "malformed Dimont model: missing bayesianNetworkSF/trees."),
    )

    # Parse each position's parameter tree
    context_positions_list = Vector{Vector{Int}}()
    nodes = Vector{DimontTreeNode}()

    for pos_elem in trees.children
        pos_elem.tag != "pos" && continue
        pt = xml_find(pos_elem, "parameterTree")
        pt === nothing &&
            throw(ModelFormatError(path, "malformed Dimont model: missing parameterTree."))

        # Parse context positions
        ctx_positions = Int[]
        ctx_poss = xml_find(pt, "contextPoss")
        if ctx_poss !== nothing
            for cp_child in ctx_poss.children
                if cp_child.tag == "pos"
                    val = _xml_numeric(cp_child, path, "Dimont contextPoss position")
                    push!(ctx_positions, Int(val))
                end
            end
        end

        # Parse the root tree element
        root_elem = xml_find(pt, "root/treeElement")
        root_elem === nothing && throw(
            ModelFormatError(path, "malformed Dimont model: missing root treeElement.")
        )

        node = _parse_tree_element(root_elem, path)
        push!(context_positions_list, ctx_positions)
        push!(nodes, node)
    end

    isempty(nodes) &&
        throw(ModelFormatError(path, "no parameter trees found in Dimont model."))

    length_val = length(nodes)
    length_val > DIMONT_MAX_LENGTH && throw(
        ModelFormatError(
            path, "Dimont length $length_val exceeds limit $DIMONT_MAX_LENGTH."
        ),
    )

    # Compute span
    span = _dimont_span(context_positions_list)
    span > DIMONT_MAX_SPAN &&
        throw(ModelFormatError(path, "Dimont span $span exceeds limit $DIMONT_MAX_SPAN."))

    # Materialize the dense tensor and flatten to 2D
    rep = _build_dimont_representation(
        nodes, context_positions_list, span, length_val, path
    )

    name = _basename_without_extension(path)
    return Dimont(name, rep, span, length_val)
end

function _parse_tree_element(elem::XMLElement, path::AbstractString)
    ctx_pos = Int(
        _xml_numeric(xml_find(elem, "contextPos"), path, "Dimont tree contextPos")
    )

    # Check for leaf (pars with pos children) vs internal (children).
    # In Jstacs XML, treeElement nodes may have both <pars> and <children> tags,
    # but <pars> only contains leaf parameters when it has <pos> children.
    pars = xml_find(elem, "pars")
    if pars !== nothing
        pars_pos_children = XMLElement[]
        for c in pars.children
            if c.tag == "pos"
                push!(pars_pos_children, c)
            end
        end
        if !isempty(pars_pos_children)
            scores = fill(-Inf, 4)
            for par_pos in pars_pos_children
                par = xml_find(par_pos, "parameter")
                par === nothing && throw(
                    ModelFormatError(path, "malformed Dimont tree: missing parameter.")
                )
                symbol = Int(
                    _xml_numeric(
                        xml_find(par, "symbol"), path, "Dimont tree parameter symbol"
                    ),
                )
                value = _xml_numeric(
                    xml_find(par, "value"), path, "Dimont tree parameter value"
                )
                if symbol < 0 || symbol > 3
                    throw(
                        ModelFormatError(
                            path, "Dimont tree parameter symbol out of range: $symbol."
                        ),
                    )
                end
                scores[symbol + 1] = value
            end
            return DimontTreeNode(ctx_pos, scores, nothing)
        end
    end

    children_elem = xml_find(elem, "children")
    if children_elem !== nothing
        children = Vector{DimontTreeNode}(undef, 4)
        fill!(children, DimontTreeNode(-1, nothing, nothing))
        for child_pos in children_elem.children
            child_pos.tag != "pos" && continue
            child_idx_str = xml_attribute(child_pos, "val")
            child_idx_str === nothing &&
                throw(ModelFormatError(path, "Dimont tree child missing val attribute."))
            child_idx = parse(Int, child_idx_str)
            if child_idx < 0 || child_idx > 3
                throw(
                    ModelFormatError(
                        path, "Dimont tree child index out of range: $child_idx."
                    ),
                )
            end
            child_elem = xml_find(child_pos, "treeElement")
            child_elem === nothing && throw(
                ModelFormatError(path, "malformed Dimont tree: missing child treeElement."),
            )
            children[child_idx + 1] = _parse_tree_element(child_elem, path)
        end
        if any(
            c -> c.context_pos == -1 && c.scores === nothing && c.children === nothing,
            children,
        )
            throw(ModelFormatError(path, "malformed Dimont tree: expected 4 children."))
        end
        return DimontTreeNode(ctx_pos, nothing, children)
    end

    return throw(
        ModelFormatError(path, "malformed Dimont tree: expected pars or children.")
    )
end

function _dimont_span(context_positions_list::Vector{Vector{Int}})
    span = 0
    for (position, positions) in enumerate(context_positions_list)
        # position is 1-indexed; Python uses 0-indexed
        pt_pos = position - 1
        if !isempty(positions)
            span = max(span, maximum(pt_pos .- positions))
        end
    end
    return span
end

function _build_dimont_representation(
    nodes::Vector{DimontTreeNode},
    context_positions_list::Vector{Vector{Int}},
    span::Int,
    length_val::Int,
    path::AbstractString,
)
    n_rows = 5^(span + 1)
    rep = Matrix{Float32}(undef, n_rows, length_val)

    # Generate all concrete ACGT contexts of length `span`
    full_contexts = _iter_contexts(span)

    for (pos_idx, node) in enumerate(nodes)
        pt_pos = pos_idx - 1  # 0-indexed position
        ctx_positions = context_positions_list[pos_idx]

        # For each concrete context, walk the tree to get scores
        # context_scores[context_key] = Vector{Float64}(4)
        context_scores = Dict{Tuple{Vararg{Int}},Vector{Float64}}()
        for ctx in full_contexts
            current = node
            while current.scores === nothing
                # Look up the parent position for this tree node
                parent_pos = current.context_pos
                ctx_val = _context_value(ctx, span, pt_pos, parent_pos, path)
                current = current.children[ctx_val + 1]
            end
            context_scores[ctx] = current.scores .+ LOG_UNIFORM_BASE
        end

        # Build the dense position tensor
        col = _build_position_column(context_scores, span, pos_idx, rep)
    end

    return rep
end

function _iter_contexts(span::Int)
    if span == 0
        return NTuple{0,Int}[()]
    end
    return vec(collect(Iterators.product(ntuple(_ -> 0:3, span)...)))
end

function _context_value(
    full_context::Tuple,
    span::Int,
    position::Int,
    absolute_position::Int,
    path::AbstractString,
)
    if absolute_position < 0
        throw(
            ModelFormatError(
                path,
                "model references position $absolute_position before the motif start at position $position.",
            ),
        )
    end
    axis = absolute_position - (position - span)
    if axis < 0 || axis >= span
        throw(
            ModelFormatError(
                path,
                "context position $absolute_position at motif position $position does not fit span $span.",
            ),
        )
    end
    return full_context[axis + 1]
end

function _build_position_column(
    context_scores::Dict{K,Vector{Float64}}, span::Int, pos_idx::Int, rep::Matrix{Float32}
) where {K}
    n_rows = 5^(span + 1)

    # Fill the column for all 5-ary codes
    for row in 1:n_rows
        code = row - 1  # 0-indexed
        digits = _decode_5ary(code, span + 1)

        # Check if any digit is N (4)
        if any(d -> d == 4, digits)
            # N-state: will be filled with minimum below
            rep[row, pos_idx] = 0.0f0  # placeholder
        else
            # Concrete ACGT: look up the score
            ctx = Tuple(digits[1:span])  # first `span` digits are context
            symbol = digits[span + 1]    # last digit is the symbol
            scores = context_scores[ctx]
            rep[row, pos_idx] = Float32(scores[symbol + 1])
        end
    end

    # Fill N-states on each context axis with the minimum over concrete nucleotides
    _fill_n_states!(rep, pos_idx, span)

    return nothing
end

function _fill_n_states!(rep::Matrix{Float32}, pos_idx::Int, span::Int)
    n_rows = 5^(span + 1)

    # Step 1: Fill context axes (0..span-1).
    # For each entry where a context axis is N (4), fill with the minimum
    # over concrete nucleotides at that axis. Only fill entries where the
    # symbol axis is also concrete (0..3); entries with symbol=N are left
    # for step 2 (matching Python's temp tensor which has no N on symbol axis).
    for axis in 0:(span - 1)
        for row in 1:n_rows
            code = row - 1
            digits = _decode_5ary(code, span + 1)
            digits[axis + 1] != 4 && continue
            digits[span + 1] == 4 && continue  # skip symbol N, fill in step 2

            min_val = Inf32
            for nt in 0:3
                probe_digits = copy(digits)
                probe_digits[axis + 1] = nt
                probe_code = _encode_5ary(probe_digits)
                val = rep[probe_code + 1, pos_idx]
                if val < min_val
                    min_val = val
                end
            end
            rep[row, pos_idx] = min_val
        end
    end

    # Step 2: Fill symbol axis (span).
    # For each entry where the symbol is N (4), fill with the minimum
    # over concrete symbols (0..3) at that position. Context axes may
    # include N-values filled in step 1, which is correct.
    for row in 1:n_rows
        code = row - 1
        digits = _decode_5ary(code, span + 1)
        digits[span + 1] != 4 && continue

        min_val = Inf32
        for nt in 0:3
            probe_digits = copy(digits)
            probe_digits[span + 1] = nt
            probe_code = _encode_5ary(probe_digits)
            val = rep[probe_code + 1, pos_idx]
            if val < min_val
                min_val = val
            end
        end
        rep[row, pos_idx] = min_val
    end
end

# _decode_5ary is defined in bamm_reader.jl and reused here.

function _encode_5ary(digits::Vector{Int})
    code = 0
    for d in digits
        code = code * 5 + d
    end
    return code
end

"""
    _xml_numeric(elem::Union{Nothing,XMLElement}, path, label)

Extract the last numeric scalar from an XML element's text content.
Returns Float64. Throws ModelFormatError if element is missing or no numeric value found.
"""
function _xml_numeric(
    elem::Union{Nothing,XMLElement}, path::AbstractString, label::AbstractString
)
    elem === nothing && throw(ModelFormatError(path, "malformed XML: missing $label."))
    text = xml_text(elem)
    # Find the last numeric value in the text
    tokens = split(text)
    for tok in reverse(tokens)
        v = tryparse(Float64, String(tok))
        v !== nothing && return v
    end
    return throw(ModelFormatError(path, "malformed XML: no numeric value for $label."))
end
