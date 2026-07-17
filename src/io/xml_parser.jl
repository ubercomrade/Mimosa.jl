# Minimal XML parser for Jstacs model files (Dimont, Slim).
#
# This is a targeted parser for the specific XML format used by Jstacs-based
# motif discovery tools. It is NOT a general-purpose XML parser — it handles
# the subset of XML needed for Dimont/Slim model files: tags with attributes,
# nested elements, text content, and comments. It does not support namespaces,
# CDATA sections, processing instructions, or external DTDs.
#
# The parser produces a lightweight DOM of XMLElement nodes that support
# ElementTree-like navigation: `xml_find`, `xml_findall`, `xml_text`,
# `xml_attribute`, and `xml_tag`.

"""
    XMLElement

A lightweight DOM node representing an XML element.

# Fields
- `tag::String`: element tag name (without namespace prefixes)
- `attributes::Vector{Pair{String,String}}`: attribute key-value pairs
- `children::Vector{XMLElement}`: child elements (in document order)
- `text::String`: text content directly inside this element (concatenated)
"""
struct XMLElement
    tag::String
    attributes::Vector{Pair{String,String}}
    children::Vector{XMLElement}
    text::String
end

"""
    xml_attribute(elem::XMLElement, key::AbstractString)

Return the attribute value for `key`, or `nothing` if not present.
"""
function xml_attribute(elem::XMLElement, key::AbstractString)
    for (k, v) in elem.attributes
        k == key && return v
    end
    return nothing
end

"""
    xml_text(elem::XMLElement)

Return the concatenated text content of `elem` (including descendant text).
"""
xml_text(elem::XMLElement) = elem.text

"""
    xml_find(elem::XMLElement, path::AbstractString)

Find the first descendant of `elem` matching the slash-separated `path`.

Path segments match direct children by tag name. A leading `//` searches
all descendants; a simple `tag1/tag2` navigates level by level.

Returns `nothing` if no match is found.
"""
function xml_find(elem::XMLElement, path::AbstractString)
    parts = split(path, '/')
    # Skip leading empty parts from // prefix
    start_idx = 1
    if isempty(parts[1])
        start_idx = 2
        # Skip any additional empty parts (from // prefix)
        while start_idx <= length(parts) && isempty(parts[start_idx])
            start_idx += 1
        end
    end

    if start_idx > 1
        # Path started with // — search all descendants
        current = _all_descendants(elem)
    else
        current = XMLElement[elem]
    end

    for i in start_idx:length(parts)
        target = parts[i]
        next_elems = XMLElement[]
        for parent in current
            for child in parent.children
                if child.tag == target
                    push!(next_elems, child)
                end
            end
        end
        isempty(next_elems) && return nothing
        current = next_elems
    end

    return isempty(current) ? nothing : current[1]
end

function _all_descendants(elem::XMLElement)
    result = XMLElement[]
    for child in elem.children
        push!(result, child)
        append!(result, _all_descendants(child))
    end
    return result
end

# ── Parser ────────────────────────────────────────────────────────────────

const XML_MAX_FILE_SIZE = 50_000_000  # 50 MB limit

"""
    _starts_at(content::String, pos::Int, prefix::AbstractString)

Check if `content` has `prefix` starting at byte position `pos`.
"""
function _starts_at(content::String, pos::Int, prefix::AbstractString)
    n = ncodeunits(prefix)
    pos + n - 1 <= ncodeunits(content) || return false
    return content[pos:(pos + n - 1)] == prefix
end

"""
    xml_parse(path::AbstractString)

Parse an XML file and return the root [`XMLElement`].

Throws `ModelFormatError` on malformed input or file size violations.
"""
function xml_parse(path::AbstractString)
    isfile(path) || throw(ModelFormatError(path, "file not found."))
    filesize = stat(path).size
    filesize > XML_MAX_FILE_SIZE && throw(
        ModelFormatError(
            path, "file size $filesize bytes exceeds limit $XML_MAX_FILE_SIZE."
        ),
    )
    content = read(path, String)
    return _xml_parse_string(content, path)
end

function _xml_parse_string(content::String, path::AbstractString)
    len = ncodeunits(content)
    pos = 1

    # Skip XML declaration and whitespace
    pos = _xml_skip_whitespace_and_decls(content, pos, len)

    # Parse the root element
    root, pos = _xml_parse_element(content, pos, len, path)
    root === nothing && throw(ModelFormatError(path, "no root element found."))
    pos = _xml_skip_whitespace_and_decls(content, pos, len)
    pos > len || throw(ModelFormatError(path, "unexpected content after root element."))
    return root
end

function _xml_skip_whitespace_and_decls(content::String, pos::Int, len::Int)
    while pos <= len
        c = content[pos]
        if isspace(c)
            pos += 1
        elseif _starts_at(content, pos, "<?xml") || _starts_at(content, pos, "<?")
            # Skip XML declaration or processing instruction
            close_idx = findnext("?>", content, pos)
            close_idx === nothing && break
            pos = last(close_idx) + 1
        elseif _starts_at(content, pos, "<!--")
            # Skip comment
            close_idx = findnext("-->", content, pos)
            close_idx === nothing && break
            pos = last(close_idx) + 1
        elseif _starts_at(content, pos, "<!DOCTYPE")
            # Skip DOCTYPE
            close_idx = findnext(">", content, pos)
            close_idx === nothing && break
            pos = last(close_idx) + 1
        else
            break
        end
    end
    return pos
end

function _xml_parse_element(content::String, pos::Int, len::Int, path::AbstractString)
    # Find opening '<'
    while pos <= len && isspace(content[pos])
        pos += 1
    end
    pos > len && return nothing, pos

    content[pos] == '<' || return nothing, pos
    pos += 1

    # Check for comment or other special tags
    if pos <= len && content[pos] == '!'
        # Skip comment or declaration
        close_idx = findnext(">", content, pos)
        close_idx === nothing && throw(ModelFormatError(path, "unclosed special tag."))
        return nothing, last(close_idx) + 1
    end

    # Parse tag name
    tag_start = pos
    while pos <= len && !_is_xml_delim(content[pos])
        pos += 1
    end
    tag = content[tag_start:(pos - 1)]
    isempty(tag) && throw(ModelFormatError(path, "empty tag name."))

    # Parse attributes
    attrs = Vector{Pair{String,String}}()
    while pos <= len
        while pos <= len && isspace(content[pos])
            pos += 1
        end
        if pos > len
            throw(ModelFormatError(path, "unclosed tag <$tag>."))
        end
        if content[pos] == '/'
            # Self-closing tag
            pos += 1
            (pos <= len && content[pos] == '>') ||
                throw(ModelFormatError(path, "malformed self-closing tag <$tag/>."))
            pos += 1
            return XMLElement(tag, attrs, XMLElement[], ""), pos
        end
        if content[pos] == '>'
            pos += 1
            break
        end
        # Parse attribute name
        attr_name_start = pos
        while pos <= len && content[pos] != '=' && !_is_xml_delim(content[pos])
            pos += 1
        end
        attr_name = content[attr_name_start:(pos - 1)]
        if pos > len || content[pos] != '='
            throw(ModelFormatError(path, "attribute without value in <$tag>."))
        end
        pos += 1
        # Skip whitespace
        while pos <= len && isspace(content[pos])
            pos += 1
        end
        if pos > len || (content[pos] != '"' && content[pos] != '\'')
            throw(ModelFormatError(path, "attribute value not quoted in <$tag>."))
        end
        quote_char = content[pos]
        pos += 1
        val_start = pos
        while pos <= len && content[pos] != quote_char
            pos += 1
        end
        if pos > len
            throw(ModelFormatError(path, "unclosed attribute value in <$tag>."))
        end
        attr_value = content[val_start:(pos - 1)]
        pos += 1
        push!(attrs, attr_name => attr_value)
    end

    # Parse content (children and text)
    children = XMLElement[]
    text_parts = String[]
    closed = false

    while pos <= len
        # Check for closing tag
        if _starts_at(content, pos, "</")
            pos += 2
            close_tag_start = pos
            while pos <= len && content[pos] != '>'
                pos += 1
            end
            if pos > len
                throw(ModelFormatError(path, "unclosed </$tag>."))
            end
            close_tag = strip(content[close_tag_start:(pos - 1)])
            close_tag != tag && throw(
                ModelFormatError(path, "mismatched tags: <$tag> closed by </$close_tag>."),
            )
            pos += 1
            closed = true
            break
        elseif _starts_at(content, pos, "<!--")
            # Skip comment
            close_idx = findnext("-->", content, pos)
            close_idx === nothing && throw(ModelFormatError(path, "unclosed comment."))
            pos = last(close_idx) + 1
        elseif content[pos] == '<'
            # Child element
            child, pos = _xml_parse_element(content, pos, len, path)
            if child !== nothing
                push!(children, child)
            end
        else
            # Text content
            text_start = pos
            while pos <= len && content[pos] != '<'
                pos += 1
            end
            text_chunk = content[text_start:(pos - 1)]
            stripped = strip(text_chunk)
            if !isempty(stripped)
                push!(text_parts, String(stripped))
            end
        end
    end

    closed || throw(ModelFormatError(path, "unclosed <$tag>."))
    return XMLElement(tag, attrs, children, join(text_parts, " ")), pos
end

function _is_xml_delim(c::Char)
    return c == ' ' ||
           c == '\t' ||
           c == '\n' ||
           c == '\r' ||
           c == '>' ||
           c == '/' ||
           c == '='
end
