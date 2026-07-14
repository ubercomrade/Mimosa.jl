# Extending Mimosa

Mimosa.jl provides a small public contract for adding custom motif models
without modifying the package. A custom model that needs no context can
participate in `scan`, `prepare_profile`, `compare`, `selectsites`, and
`reconstruct_pfm` by implementing only three methods. Models with
context implement two additional geometry methods.

The contract is defined by ADR 0003 (`docs/adr/0003-model-geometry-contract.md`
at the repository root) and the Extensibility API Plan
(`Mimosa.jl/EXTENSIBILITY_API_PLAN.md`).

## Minimal contract for comparison

A custom model subtypes `Mimosa.AbstractMotifModel` and implements:

```julia
import Mimosa

struct MyModel <: Mimosa.AbstractMotifModel
    label::String          # any field names you like; no `name`/`weights`/
                          # `representation`/`order`/`span` requirement
    pattern::Vector{UInt8} # encoded consensus (A=0,C=1,G=2,T=3,N=4)
end

Mimosa.modelname(model::MyModel) = model.label
Mimosa.motif_length(model::MyModel) = length(model.pattern)

function Mimosa.scan_pair_kernel!(
    fwd_out::AbstractVector{Float32},
    rev_out::AbstractVector{Float32},
    model::MyModel,
    sequence::AbstractVector{UInt8},
    n_positions::Int,
)
    pat = model.pattern
    rc_pat = UInt8[b == 0x04 ? b : (0x03 - b) for b in reverse(pat)]
    @inbounds for pos in 1:n_positions
        f = zero(Float32)
        r = zero(Float32)
        for k in 1:length(pat)
            b = sequence[pos + k - 1]
            f += (b == pat[k]) ? 1.0f0 : 0.0f0
            r += (b == rc_pat[k]) ? 1.0f0 : 0.0f0
        end
        fwd_out[pos] = f
        rev_out[pos] = r
    end
    return (fwd_out, rev_out)
end
```

That is the entire contract for `:compare`. After these definitions:

```julia
sequences = Mimosa.readsequences("sequences.fasta")[1]
result = compare(MyModel(...), MyModel(...), sequences; metric=:co)
```

works through the same scan → normalize → anchor → profile alignment
path as the built-in models.

## Geometry: models with context

If the model needs `L` bases before the motif site or `R` bases after it
to compute one score, override the corresponding context accessor:

```julia
Mimosa.left_context(model::MyModel) = L
Mimosa.right_context(model::MyModel) = R
```

Both default to `0`. The derived quantities are computed by Mimosa.jl:

| Quantity | Definition |
|---|---|
| `window_size(model)` | `left_context(model) + motif_length(model) + right_context(model)` |
| `npositions(model, sequence_length)` | `max(sequence_length - window_size(model) + 1, 0)` |
| `site_start_offset(model)` | `left_context(model)` |

Forward and reverse scores at the same scan index share the same
physical window and the same physical site interval. The reverse
kernel is responsible for orienting the score computation. If a future
model needs different physical site intervals per orientation, that is
outside this contract and requires a separate ADR.

## Interface validation

`Mimosa.validate_model(model; capability=:compare)` returns the model
on success and throws `Mimosa.ModelInterfaceError` on failure. Supported
capabilities:

| Capability | Required methods |
|---|---|
| `:compare` | `AbstractMotifModel` subtype, `modelname`, `motif_length`, `scan_pair_kernel!` |
| `:sites` | `:compare` plus valid geometry (positive `motif_length`, non-negative contexts) |
| `:cache` | `:compare` plus `model_fingerprint` returning a non-empty SHA-256 string |

The public `scan`, `prepare_profile`, `compare`, `selectsites`, and
cache/null entry points perform the applicable validation at the boundary;
inner loops and worker tasks do not repeat it.

## Optional capabilities

- `scorebounds(model)` is needed only by `inspect-model` and by formats
  that export theoretical score bounds. Plain `compare` does not call
  it.
- `model_fingerprint(model)` is needed only by the cache and null
  bundle compatibility tracking. Plain `compare` does not call it. A
  custom model implements it as a stable SHA-256 hex string of all
  parameters that affect scores.
- Specialized `scan_forward!`/`scan_reverse!`/`scan_best!`/`scan_both!`
  methods are optional performance overrides; the generic fallback
  already drives `scan_pair_kernel!` correctly. Add one only after a
  benchmark justifies it.

## External score adapter (no `AbstractMotifModel`)

If you have an external scanner and do not need site extraction or the
built-in batch scanning, use `ScoreProfile`:

```julia
forward_scores = Mimosa.scan(my_model, batch; strands=ForwardOnly())  # or external scores
profile = Mimosa.ScoreProfile("name", forward_scores)
prepared = Mimosa.prepare_profile(profile)
result = Mimosa.compare(prepared, other_prepared; metric=:co)
```

`ScoreProfile(name, scores)` constructs a symmetric profile where both
strands use the same scores.

## File formats

Mimosa.jl keeps file parsing separate from model algorithms. Built-in readers
(`readmodel`, `read_bamm`, `read_sitega`, `read_dimont`, `read_slim`,
`read_scores`) cover the supported text and bundle formats. The
keyword `format=:auto` remains a symbol-based boundary adapter for
built-in formats. A typed third-party format extension API is planned
for Stage 6 but is not part of the current public contract. External
packages should parse their formats themselves and construct a custom
model or `ScoreProfile` through the public constructors.

The portable bundle format stores only built-in model kinds. A future
custom-bundle codec protocol is a separate versioned decision and is
not part of this extension contract.

## Principles

- **No central registry**: models and formats are added via dispatch,
  not by editing a registry dictionary.
- **No structural duck typing**: generic workflows call
  `modelname(model)`, `motif_length(model)`, `left_context(model)`,
  `right_context(model)`, and `scan_pair_kernel!`. They never access
  `model.name`, `model.representation`, `model.weights`, `model.order`,
  or `model.span` on an `AbstractMotifModel` value.
- **No `Any` fields**: all model fields must have concrete or
  parametric types.
- **No string dispatch in hot paths**: strings and symbols are parsed
  only at API or I/O boundaries.
- **No type piracy**: only define methods for types you own.
- **Float32 scan values, deterministic tie-breaking, and one-based
  inclusive coordinates are preserved** by every extension.

## Testing custom models

A downstream test module should:

- define the model in a separate module that imports Mimosa as a
  regular dependency;
- not access `Mimosa._private_name`;
- not have fields named `representation`, `weights`, `order`, or `span`;
- cover `validate_model` success and failure for each required method;
- cover single/batch/empty/too-short/exact-window sequences;
- cover all four strand policies;
- cover serial/threaded exact equivalence and order preservation;
- cover worker exception propagation;
- cover scalar, prepared, and one-to-many comparison;
- cover custom/built-in comparison in both orders;
- cover sites and PFM reconstruction (with and without context);
- cover the absence of a fingerprint requirement for plain `compare`
  and the explicit error for `:cache`.

See `test/downstream/runtests.jl` for a working example.
