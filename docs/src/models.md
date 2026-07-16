# Supported Models

Mimosa.jl supports six model families, each with a concrete immutable type
and specific file format.

## Model types

| Type | Family | Format | Description |
|------|--------|--------|-------------|
| `PWM{T,M,B}` | Matrix | MEME, PFM | Position Weight Matrix with N-state row |
| `BaMM{T,M}` | Higher-order | `.ihbcp` | Bayesian Markov Model |
| `SiteGA{T,M}` | Higher-order | `.mat` | Dinucleotide model |
| `Dimont{T,M}` | Higher-order | XML | Jstacs Bayesian network |
| `Slim{T,M}` | Higher-order | XML | Jstacs GenDisMix classifier |

## Precomputed profiles

`ScoreProfile` is a precomputed profile source, not a motif model. Its `length`
is the number of profile rows for compatibility; use `nrows(profile.scores)` for
that value. It is prepared directly for profile comparison and cannot be scanned
or serialized as a motif bundle.

## File formats

### MEME (PWM)

MEME letter-probability matrices are read as frequencies and converted to a
5-row log-odds PWM (A, C, G, T, N), with the N row as the per-column minimum.

### PFM (Position Frequency Matrix) input

A plain-text `.pfm` file is an input format, not a public model type.
`readmodel("motif.pfm"; background=0.25f0)` reads its 4-row (A, C, G, T)
frequency matrix and returns a ready-to-scan `PWM`. For an in-memory frequency
matrix, call `pwm_from_pfm(pfm; background=...)`; the conversion applies the
documented pseudocount and adds the N-state row.

### BaMM `.ihbcp`

Bayesian Markov Model with higher-order context. The `.ihbcp` format stores
a flattened matrix of shape `(5^(order+1), motif_length)`.

### SiteGA `.mat`

Dinucleotide model with 25-row representation (5×5 dinucleotide codes),
flattened from Python's `(5, 5, length)` to `(25, length)`.

### Dimont XML

Jstacs XML format. The parser extracts `MarkovModelDiffSM` parameters,
materializes them into a dense 5-ary tensor, and flattens to
`(5^(span+1), motif_length)`.

### Slim XML

Jstacs XML format. The parser extracts `SLIM` component/ancestor parameters,
normalizes via log-sum-exp, and materializes into the same representation
as Dimont.

## Score bounds

Each model type implements `scorebounds(model)` returning `(min_score, max_score)`
—the per-column minimum and maximum, summed across positions. This is used
for score normalization and null distribution fitting.

## Scanning

All models support the same scanning API via multiple dispatch:

```julia
scan(model, sequence; strands=BestStrand())
scan(model, batch; strands=BestStrand(), execution=SerialExecution())
scan!(dest, model, sequence; strands=ForwardOnly())
```

The scanning kernel is shared across model families and parameterized by
the [public geometry contract](extending.md): `motif_length`, `left_context`,
and `right_context`. Mimosa.jl derives
`window_size`, `npositions`, and `site_start_offset` from these.
For BaMM, Dimont, and Slim, the window includes preceding context
(`left_context = order/span`); for PWM and SiteGA it equals the motif
length. `npositions(model, sequence_length)` exposes the resulting
scan-track length.

## Custom models

See [Extending Mimosa](extending.md). A custom model subtypes
`AbstractMotifModel` and implements `modelname`, `motif_length`, and
`scan_kernel!`. Context models additionally implement
`left_context` and/or `right_context`. The generic scan, prepare,
compare, sites, and PFM reconstruction workflows then work through the
public API without modifying Mimosa.jl.

## Type parameters

All matrix model types are parametric:

```julia
PWM{T<:AbstractFloat, M<:AbstractMatrix{T}, B<:NTuple{4,AbstractFloat}}
BaMM{T<:AbstractFloat, M<:AbstractMatrix{T}}
```

This allows the compiler to specialize on the concrete element type (typically
`Float32`) and array type (typically `Matrix{Float32}`), ensuring type-stable
hot kernels with zero per-position allocations.
