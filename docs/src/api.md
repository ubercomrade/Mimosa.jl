# Julia API Guide

This page is the practical entry point to Mimosa.jl from Julia. The usual
workflow is:

1. read motif models and DNA sequences;
2. scan models or compare them through their score profiles;
3. inspect sites, reconstruct a PFM, or estimate significance when needed.

Mimosa compares models by their behavior on the same sequences. It does not
convert BaMM, SiteGA, Dimont, or Slim parameters into a PWM before comparison.

## Install and load

Mimosa.jl requires Julia 1.12 or newer. Until it is registered in the General
registry, install it from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/ubercomrade/Mimosa.jl.git")

using Mimosa
```

## Read models and sequences

`readmodel` detects the input format from the filename. `readsequences` reads
FASTA input and returns both encoded sequences and their names:

```julia
model = readmodel("examples/pif4.meme")
sequences, names = readsequences("examples/foreground.fa")

model_name = modelname(model)
first_sequence = sequence(sequences, 1)
```

Supported model inputs are:

| Input | Reader | Result |
|---|---|---|
| MEME PWM or `.pfm` | `readmodel` | `PWM` |
| BaMM `.ihbcp` | `readmodel` or `read_bamm` | `BaMM` |
| SiteGA `.mat` | `readmodel` or `read_sitega` | `SiteGA` |
| Dimont XML | `readmodel` or `read_dimont` | `Dimont` |
| Slim XML | `readmodel` or `read_slim` | `Slim` |
| FASTA-like numeric rows | `read_scores` | `ScoreProfile` |
| Mimosa model bundle | `readmodel` | built-in model type |

FASTA bases are case-insensitive. `N` and IUPAC ambiguous bases are encoded as
`N_CODE`; they are accepted by the scanner but are not counted when a PFM is
reconstructed. See [Supported Models](models.md) for format details.

## Scan a model

Use a strand policy explicitly when the result must be reproducible and easy
to interpret:

```julia
forward = scan(model, sequences; strands=ForwardOnly())
reverse = scan(model, sequences; strands=ReverseOnly())
best = scan(model, sequences; strands=BestStrand())
both = scan(model, sequences; strands=BothStrands())

forward_first = row(forward, 1)
reverse_first = row(both.reverse, 1)
```

For a batch, `scan` returns a `RaggedArray{Float32}` because sequences can have
different lengths. `BothStrands()` returns a `StrandPair` with separate
`forward` and `reverse` ragged arrays. A single-sequence scan returns a
`Vector{Float32}` (or a `StrandPair` of vectors).

The scan track contains one value per complete model window. For PWM and
SiteGA, a window is the motif itself. BaMM, Dimont, and Slim use context on
both sides of the motif. The returned positions are one-based scan positions;
for a context model, the motif starts at
`scan_position + site_start_offset(model)`. See [Data Layout](data_layout.md)
for coordinate details.

`SerialExecution()` is the default. Threaded scanning requires Julia runtime
threads and an explicit policy:

```bash
JULIA_NUM_THREADS=4 julia --project=.
```

```julia
scores = scan(
    model,
    sequences;
    strands=BestStrand(),
    execution=ThreadedExecution(4),
)
```

Serial and threaded execution preserve result order and exact discrete fields.

For profile comparisons and `build_null`, `outer_execution` is the outer level
(targets or null pairs) and `scan_execution` controls scanning sequences within
one profile. Choose one threaded level: using multi-threaded policies for both
throws `ArgumentError`. For one scalar comparison use
`scan_execution=ThreadedExecution(n)`; for many targets or null pairs use
`outer_execution=ThreadedExecution(n)` with the default
serial `scan_execution`.

## Compare motif models

Comparison requires the same `EncodedSequenceBatch` for both models:

```julia
query = readmodel("examples/pif4.meme")
target = readmodel("examples/gata2.meme")
sequences, _ = readsequences("examples/foreground.fa")

result = compare(
    query,
    target,
    sequences;
    metric=:co,
    search_range=10,
    window_radius=10,
    realign_window=3,
)

result.score       # higher means a stronger match
result.offset      # query displacement relative to target
result.orientation # "++", "+-", "-+", or "--"
result.n_sites
```

The pipeline scans both models, normalizes scores to empirical `-log10(FPR)`
values, selects anchors, and compares windows across shifts and orientations.
A positive `offset` means that the query is shifted to the right.

Hybrid empirical normalization with an exact high-score tail is the standard
normalization used by the Julia API. The exact full-table implementation is an
internal reference implementation rather than a selectable public strategy.

Available metrics are:

| Metric | Meaning |
|---|---|
| `:co` | pooled continuous overlap |
| `:co_rowwise` | overlap averaged so each anchor contributes equally |
| `:dice` | pooled Dice similarity |
| `:dice_rowwise` | Dice averaged per anchor |
| `:cosine` | cosine similarity averaged per anchor |

Use a separate FASTA batch for normalization when the comparison sequences are
not an appropriate background:

```julia
background, _ = readsequences("examples/background.fa")
result = compare(query, target, sequences; background, metric=:cosine)
```

`ComparisonResult` is the Julia result object. JSON serialization is available
for comparison results and annotated results, not for `SiteCollection`:

```julia
println(to_json(result))
payload = to_dict(result)
```

## Reuse prepared profiles

When one query is compared with many targets, prepare it once. Preparation
performs scanning, normalization, and anchor collection:

```julia
cache = Cache(".mimosa-cache")
prepared_query = prepare_profile(query, sequences; min_logfpr=0.0, cache=cache)
results = compare(
    prepared_query,
    [target, readmodel("examples/foxa2.meme")],
    sequences;
    metric=:cosine,
    outer_execution=ThreadedExecution(4),
    cache=cache,
)
```

`cache` is opt-in and stores normalized profiles plus anchors. Cache keys include
the model content, comparison sequences, optional normalization background, and
`min_logfpr`; changing any of these inputs produces a cache miss. Use
`mimosa cache clear --cache-dir .mimosa-cache` to remove entries.

For already computed score rows, use `ScoreProfile` and `read_scores`:

```julia
query_scores = read_scores("examples/scores_1.fasta")
target_scores = read_scores("examples/scores_2.fasta")
prepared_scores = prepare_profile(query_scores)
score_result = compare(prepared_scores, target_scores; metric=:co)
```

Prepared profiles compared directly must use the same `min_logfpr` threshold.

## Extract sites and reconstruct a PFM

Selectors determine which hits are retained:

```julia
best_sites = selectsites(
    model,
    sequences,
    BestPerSequence();
    strands=BestStrand(),
)

threshold_sites = selectsites(
    model,
    sequences,
    ThresholdHits(5.0f0);
    strands=BothStrands(),
)

top_sites = selectsites(model, sequences, TopFractionHits(0.1))
```

`SiteCollection` stores parallel arrays: `seq_indices`, `starts`, `strands`,
and `scores`. `starts` are one-based scan positions, and `strands` are `0`
(forward) or `1` (reverse). For models with context, add
`site_start_offset(model)` to a start to get the physical motif start.

Extract canonical-orientation site strings when needed:

```julia
site_matrix = extract_site_matrix(
    sequences,
    best_sites,
    motif_length(model);
    site_offset=site_start_offset(model),
)
strings = site_strings(site_matrix)
```

Reconstruct a 4-row (A, C, G, T) PFM directly:

```julia
pfm = reconstruct_pfm(
    model,
    sequences,
    BestPerSequence();
    strands=BestStrand(),
    pseudocount=1.0f-4,
)
```

`N` bases are skipped in the counts. The reverse strand is reverse-complemented
before counting, so PFM columns are in canonical motif orientation.

## Significance with a null distribution

A profile null samples distinct ordered model pairs with replacement and fits a
GEV distribution to their comparison scores. PWM models can be independently
shuffled before every comparison:

```julia
models = [query, target, readmodel("examples/foxa2.meme")]

built = build_null(
    models,
    sequences,
    metric=:co,
    n_samples=2000,
    shuffle=true,
    seed=127,
)
savenull("output/null_bundle", built.distribution)

distribution = loadnull("output/null_bundle")
annotated = annotate_results([result], distribution)
println(to_json(first(annotated)))
```

`annotate_results` adds p-values, Benjamini–Hochberg adjusted p-values, and
E-values. Null bundles are profile-only and use format version 5. They are
compatible only with the same metric, comparison settings, sequences,
background, and model family. The manifest records the source-model fingerprint,
sampling seed, shuffle flag, and sampling algorithm version. See
[Storage Format](storage.md) for the portable layout.

## Errors and troubleshooting

Input and model failures use typed errors:

```julia
try
    model = readmodel("input.meme")
catch error
    error isa MimosaError || rethrow()
    println("Mimosa input error: ", error)
end
```

Common fixes:

- **No scan positions**: the sequence is shorter than the model window.
- **Different prepared thresholds**: prepare all profiles with the same
  `min_logfpr`.
- **No null targets**: check motif names and ensure each query has an eligible
  motif from another group.
- **Too few Julia threads**: set `JULIA_NUM_THREADS` before starting Julia;
  `ThreadedExecution` cannot create runtime threads.

## Core API entries

The following entries document the objects used in the workflows above. The
page intentionally starts with tasks and examples rather than an alphabetical
export list.

```@docs
readmodel
read_scores
EncodedSequenceBatch
ScoreProfile
scan
StrandPolicy
ForwardOnly
ReverseOnly
BestStrand
BothStrands
StrandPair
compare
ComparisonResult
prepare_profile
PreparedProfile
selectsites
SiteCollection
reconstruct_pfm
NullDistribution
NullBuildResult
AnnotatedResult
annotate_results
to_dict
to_json
ExecutionPolicy
SerialExecution
ThreadedExecution
MimosaError
```
