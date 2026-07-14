# Small, meaningful exception hierarchy for Mimosa.

"""
    MimosaError

Abstract supertype of all Mimosa-specific exceptions.
"""
abstract type MimosaError <: Exception end

"""
    ModelFormatError

Raised when a model file is malformed, has an unsupported format, or contains
invalid data. `path` is the file path (may be empty for in-memory errors).
"""
struct ModelFormatError <: MimosaError
    path::String
    message::String
end

"""
    ModelDimensionError

Raised when a model has incompatible or invalid dimensions.
"""
struct ModelDimensionError <: MimosaError
    message::String
end

"""
    InvariantError

Raised when an internal invariant is violated.
"""
struct InvariantError <: MimosaError
    message::String
end

"""
    ModelInterfaceError

Raised when a model fails the public extension interface validation
performed by [`validate_model`](@ref) or by the public `scan`,
`prepare_profile`, `compare`, `selectsites`, and cache/null entry points.

Fields:
- `capability::Symbol`: capability being validated (`:compare`, `:sites`,
  `:cache`).
- `model_type::String`: `typeof(model)` string for diagnostics.
- `message::String`: human-readable description of the violation.
"""
struct ModelInterfaceError <: MimosaError
    capability::Symbol
    model_type::String
    message::String
end

function Base.showerror(io::IO, e::ModelFormatError)
    return print(io, "ModelFormatError: $(e.path): $(e.message)")
end

function Base.showerror(io::IO, e::ModelDimensionError)
    return print(io, "ModelDimensionError: $(e.message)")
end

function Base.showerror(io::IO, e::InvariantError)
    return print(io, "InvariantError: $(e.message)")
end

function Base.showerror(io::IO, e::ModelInterfaceError)
    return print(io, "ModelInterfaceError ($(e.capability), $(e.model_type)): $(e.message)")
end
