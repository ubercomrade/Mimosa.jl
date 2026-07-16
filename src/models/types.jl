# Abstract type hierarchy for Mimosa models and profile sources.
#
# Defined before geometry.jl so the public geometry accessors can refer to
# `AbstractMotifModel` in their default method signatures.

"""
    AbstractProfileSource

Abstract supertype of inputs that can be prepared for profile comparison.
"""
abstract type AbstractProfileSource end

"""
    AbstractMotifModel

Abstract supertype of all motif model families (PWM, PFM, BaMM, SiteGA, etc.).
Motif models can be scanned against encoded sequences. Precomputed profiles are
`AbstractProfileSource`s, but are not motif models.

A custom model subtypes `AbstractMotifModel` and implements the minimal
contract described in `docs/src/extending.md`:

- `modelname(model)::AbstractString`
- `motif_length(model)::Integer`
- `scan_kernel!(forward, reverse, model, sequence, n_positions)` (or
  specialized `scan_forward!`/`scan_reverse!`/`scan_best_strand!`/`scan_both!`
  methods)

`left_context` and `right_context` default to zero. `window_size`,
`npositions`, and `site_start_offset` are derived by Mimosa.jl.
"""
abstract type AbstractMotifModel <: AbstractProfileSource end

abstract type AbstractContextModel{T<:AbstractFloat} <: AbstractMotifModel end

"""Return whether `model` has a direct sequence-scanning implementation.

This is a transitional compatibility shim. Subtyping
`AbstractMotifModel` is itself the scannability declaration; new model
types do not need to define `is_scannable` and the function may be
deprecated in a future breaking release.
"""
is_scannable(::AbstractMotifModel) = true
