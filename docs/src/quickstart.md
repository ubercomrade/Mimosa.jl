# Quick Start

## Installation

Mimosa.jl requires Julia 1.12 or newer. Until it is registered in General,
install it from its repository:

```julia
using Pkg
Pkg.add(url="https://github.com/ubercomrade/Mimosa.jl.git")
```

After registration, use `Pkg.add("Mimosa")`. The examples below assume the
repository root is the working directory.

## Read and Scan a Model

```julia
using Mimosa

model = readmodel("examples/pif4.meme")
sequences, names = readsequences("examples/foreground.fa")

scores = scan(model, sequences; strands=BestStrand())
```

`EncodedSequenceBatch` and `RaggedArray` use flat offset-based storage. Scan
policies are `ForwardOnly()`, `ReverseOnly()`, `BestStrand()`, and
`BothStrands()`.

## Compare Models

Model comparison requires sequences and follows scan, normalization, anchor,
and profile alignment:

```julia
query = readmodel("examples/pif4.meme")
target = readmodel("examples/gata2.meme")
sequences, _ = readsequences("examples/foreground.fa")

comparison = compare(
    query,
    target,
    sequences;
    metric=:co,
    search_range=10,
    window_radius=10,
    realign_window=3,
)

println(to_json(comparison))
```

Available metrics are `:co`, `:co_rowwise`, `:dice`, `:dice_rowwise`, and
`:cosine`.

## Reuse Prepared Profiles

Prepare profiles once when the same query is compared repeatedly:

```julia
query_scores = ScoreProfile("query", scores)
target_scores = ScoreProfile("target", scores)

prepared_query = prepare_profile(query_scores; min_logfpr=0.0f0)
prepared_target = prepare_profile(target_scores; min_logfpr=0.0f0)

prepared_comparison = compare(prepared_query, prepared_target; metric=:cosine)
prepared_results = compare(prepared_query, [prepared_target]; metric=:cosine)
```

Prepared profiles compared together must use the same `min_logfpr` threshold.

## Sites and PFM Reconstruction

```julia
sites = selectsites(model, sequences, BestPerSequence(); strands=BestStrand())

pfm = reconstruct_pfm(
    model,
    sequences,
    TopFractionHits(0.1);
    strands=BestStrand(),
    pseudocount=1.0f-4,
)
```

Site collections store one-based scan positions in `starts`. For models with
context, the physical motif starts at `start + site_start_offset(model)`. JSON
serialization is available for comparison results, not site collections.

## Null Distributions

```julia
models = [query, target, readmodel("examples/foxa2.meme")]

null_result = build_null(
    models;
    sequences=sequences,
    metric=:co,
    n_samples=2000,
    shuffle=true,
    seed=127,
)

dist = null_result.distribution
savenull("output/null_bundle", dist)
loaded = loadnull("output/null_bundle")
annotated = annotate_results([comparison], loaded; effective_number_of_targets=1)
```

Null bundles use format version 5 and strategy `"profile"` only. Older bundles must be rebuilt.

## Threaded Execution

Threading is opt-in in both the runtime and the API:

```bash
JULIA_NUM_THREADS=4 julia --project=.
```

Mimosa keeps all comparison work in that single Julia process. For a scalar
motif comparison, thread sequence scanning explicitly:

```julia
comparison = compare(
    query, target, sequences;
    execution=Execution(4),
)
```

Use the same `execution` setting for `build_null` and one-to-many comparisons.
Targets remain serial while each comparison uses threaded computational
kernels.

```julia
scores = scan(
    model,
    sequences;
    strands=BestStrand(),
    execution=Execution(4),
)
```

`Execution()` is the default and is equivalent to `Execution(1)`. Values
greater than one enable kernel-level parallelism, capped by the number of
Julia runtime threads. Sequential and parallel workflows preserve result
ordering and exact discrete fields.
