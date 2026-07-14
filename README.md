# Mimosa.jl

Mimosa.jl is a Julia 1.10+ package and command-line tool for DNA motif comparison. It supports PWM/PFM, BaMM, SiteGA, Dimont, Slim, and precomputed score profiles. Mimosa compares heterogeneous motif models through their behavior on the same DNA sequences rather than by forcing their internal representations into a common matrix.

## Background

Transcription factors regulate gene expression by recognizing short DNA segments called transcription factor binding sites (TFBSs). Sites bound by the same factor are similar but not identical and are therefore represented by motif models. PWMs assume independent motif positions, whereas BaMM, SiteGA, Dimont, and Slim models can describe dependencies between nucleotides.

Directly aligning the parameters of different model families is generally not meaningful. Mimosa instead uses the following profile comparison pipeline:

1. Scan sequences by motifs to get score profile.
2. Convert raw scores to empirical `-log10(FPR)` values, optionally using a
   separate background sequence set for calibration.
3. Select one best anchor per sequence when `-log10(FPR) == 0`, or all anchors
   above the requested threshold.
4. Compare site-centered windows over shifts and the four strand orientations
   `++`, `+-`, `-+`, and `--`.
5. Return the highest-scoring alignment, its offset, orientation, and number of
   contributing sites.

The reported offset is the query displacement relative to the target: a
positive value means that the query is shifted to the right. Ties are resolved
deterministically by score, contributing site count, smaller absolute shift,
and orientation order `++`, `+-`, `-+`, `--`.

### Similarity metrics

For two aligned, non-negative profile windows `v1` and `v2`, the pooled
continuous overlap and Dice metrics are

```math
CO(v_1,v_2) =
\frac{\sum_i \min(v_{1i},v_{2i})}
     {\min(\sum_i v_{1i},\sum_i v_{2i})}
```

```math
Dice(v_1,v_2) =
\frac{2\sum_i \min(v_{1i},v_{2i})}
     {\sum_i v_{1i}+\sum_i v_{2i}}.
```

`co_rowwise` and `dice_rowwise` calculate the corresponding value for each
selected window and average finite values, so each anchor contributes equally.
`cosine` is also calculated per window and averaged:

```math
Cosine(v_1,v_2) =
\frac{v_1 \cdot v_2}{\lVert v_1 \rVert_2\lVert v_2 \rVert_2}.
```

All metrics are oriented so that a larger value means a stronger match.

### Statistical significance

`build_null` compares eligible cross-group model pairs and pools their scores
into an empirical null sample. Mimosa fits a generalized extreme value (GEV)
distribution in `Float64`; the upper-tail survival probability of an observed
score is its p-value. `annotate_results` additionally computes Benjamini-Hochberg
adjusted p-values and E-values.

Null bundles are reusable only when their profile metric, comparison settings,
sequence/background fingerprints, and format contract are compatible with the
observed comparison.

## Supported inputs

| Family | CLI type | Input |
|---|---|---|
| Precomputed profiles | `scores` | FASTA-like numeric score rows |
| PWM/PFM | `pwm` | MEME, plain PFM, or a Mimosa model bundle |
| BaMM | `bamm` | `.ihbcp` or a Mimosa model bundle |
| SiteGA | `sitega` | `.mat` or a Mimosa model bundle |
| Dimont | `dimont` | Jstacs XML or a Mimosa model bundle |
| Slim | `slim` | Jstacs XML or a Mimosa model bundle |

A plain PFM is converted to a scannable PWM. A precomputed `ScoreProfile` can
be compared to another `ScoreProfile`, but mixed score-profile/motif comparison
must first prepare both inputs as profiles.

## Installation

Mimosa.jl requires Julia 1.10 or newer and is not currently registered in the
General registry. Clone the repository and develop the package from its
`Mimosa.jl/` subdirectory:

```bash
git clone https://github.com/ubercomrade/mimosa.git
cd mimosa
julia --project=Mimosa.jl -e 'using Pkg; Pkg.instantiate()'
```

To make the package available in another Julia environment:

```julia
using Pkg
Pkg.develop(path="/path/to/mimosa/Mimosa.jl")
```

Verify the installation with:

```bash
julia --project=Mimosa.jl -e 'using Mimosa; println(Base.pkgversion(Mimosa))'
```

The examples below are run from the repository root.

## Quick start

```julia
using Mimosa

query = readmodel("examples/pif4.meme")
target = readmodel("examples/gata2.meme")
sequences, names = readsequences("examples/foreground.fa")

result = compare(query, target, sequences; metric=:co)
println(to_json(result))

scores = scan(query, sequences; strands=BestStrand())
sites = selectsites(query, sequences, BestPerSequence())
pfm = reconstruct_pfm(query, sequences, TopFractionHits(0.1))
```

Julia coordinates returned by the API are one-based and inclusive. JSON output
uses zero-based, half-open site coordinates at the serialization boundary.

## Command-line interface

Run the CLI directly from the checkout:

```bash
julia --project=Mimosa.jl Mimosa.jl/app/mimosa.jl --help
```

Successful commands write JSON only to `stdout`; diagnostics and errors go to
`stderr`. Exit codes are `0` for success, `1` for invalid usage, and `2` for
input, runtime, or scientific errors.

### Compare two models

```bash
JULIA_NUM_THREADS=4 julia --project=Mimosa.jl Mimosa.jl/app/mimosa.jl \
  profile examples/pif4.meme examples/gata2.meme \
  --model1-type pwm --model2-type pwm \
  --fasta examples/foreground.fa \
  --background examples/background.fa \
  --metric co --min-logfpr 2 --threads 4
```

The available metrics are `co`, `co_rowwise`, `dice`, `dice_rowwise`, and
`cosine`. Important profile options are:

| Option | Meaning | Default |
|---|---|---|
| `--fasta` | FASTA sequences scanned by motif inputs | generated sequences |
| `--background` | separate FASTA for empirical normalization | comparison FASTA |
| `--num-sequences` | generated sequence count when FASTA is absent | `1000` |
| `--seq-length` | generated sequence length | `200` |
| `--seed` | random-sequence seed | `127` |
| `--search-range` | largest tested displacement | `10` |
| `--window-radius` | site-centered profile half-window | `10` |
| `--realign-window` | local target-anchor search radius | `3` |
| `--min-logfpr` | anchor threshold; `0` selects the best per sequence | `0` |
| `--threads` | explicit API worker budget | `1` |

Precomputed profiles do not require FASTA input:

```bash
julia --project=Mimosa.jl Mimosa.jl/app/mimosa.jl \
  profile examples/scores_1.fasta examples/scores_2.fasta \
  --model1-type scores --model2-type scores --metric cosine
```

A typical result is:

```json
{
  "query": "pif4",
  "target": "gata2",
  "score": 0.73,
  "offset": -1,
  "orientation": "+-",
  "metric": "co",
  "n_sites": 42
}
```

### Build and use a null distribution

The group table is TSV or CSV and must contain motif-name and group columns.
Pairs from different groups form the null sample.

```bash
julia --project=Mimosa.jl Mimosa.jl/app/mimosa.jl \
  build-null motifs/ \
  --model-type pwm --groups groups.tsv \
  --name-column motif --group-column group \
  --fasta examples/foreground.fa --metric co \
  --output output/null_bundle

julia --project=Mimosa.jl Mimosa.jl/app/mimosa.jl \
  profile examples/pif4.meme examples/gata2.meme \
  --model1-type pwm --model2-type pwm \
  --fasta examples/foreground.fa --metric co \
  --pvalue --null-distribution output/null_bundle
```

`--strict` makes an insufficient number of unrelated targets an error;
`--min-null-targets` controls the minimum. `--jobs` remains a deprecated alias
for `--threads` on `build-null` only.

### Inspect and convert models

```bash
julia --project=Mimosa.jl Mimosa.jl/app/mimosa.jl \
  inspect-model examples/foxa2.ihbcp --type bamm

julia --project=Mimosa.jl Mimosa.jl/app/mimosa.jl \
  convert-model examples/pif4.meme output/pif4_bundle --type pwm
```

`convert-model` writes a portable Mimosa bundle. `cache clear` removes derived
profile artifacts from a selected cache directory:

```bash
julia --project=Mimosa.jl Mimosa.jl/app/mimosa.jl \
  cache clear --cache-dir .mimosa-cache
```

## Julia API

### Read, write, and scan models

```julia
using Mimosa

pwm = readmodel("examples/pif4.meme")
bamm = readmodel("examples/foxa2.ihbcp")
sequences, names = readsequences("examples/foreground.fa")

forward = scan(pwm, sequences; strands=ForwardOnly())
best = scan(pwm, sequences; strands=BestStrand())
both = scan(pwm, sequences; strands=BothStrands())

writemodel("output/pif4_bundle", pwm)
restored = readmodel("output/pif4_bundle")
```

`EncodedSequenceBatch` stores one validated flat `UInt8` buffer plus offsets;
DNA codes are `A=0x00`, `C=0x01`, `G=0x02`, `T=0x03`, and ambiguous/N=`0x04`.
Batch scan results use flat `RaggedArray` storage and preserve empty sequences
and input order.

### Compare models and profiles

```julia
background, _ = readsequences("examples/background.fa")

result = compare(
    pwm,
    bamm,
    sequences;
    background=background,
    metric=:dice_rowwise,
    search_range=10,
    window_radius=10,
    realign_window=3,
    min_logfpr=2.0,
)
```

For repeated comparisons, prepare each model once. Prepared profiles own their
normalization and anchors, and both sides must use the same `min_logfpr`:

```julia
query = prepare_profile(pwm, sequences; background=background, min_logfpr=2.0)
targets = [
    prepare_profile(readmodel("examples/gata2.meme"), sequences;
                    background=background, min_logfpr=2.0),
    prepare_profile(readmodel("examples/foxa2.meme"), sequences;
                    background=background, min_logfpr=2.0),
]

results = compare(query, targets; metric=:co)
```

External scanners can provide a `ScoreProfile` instead of implementing a motif
model:

```julia
profile = ScoreProfile("external", external_scores) # RaggedArray{Float32}
prepared = prepare_profile(profile; min_logfpr=0.0)
result = compare(prepared, other_prepared; metric=:cosine)
```

### Extract sites and reconstruct a PFM

```julia
best_sites = selectsites(pwm, sequences, BestPerSequence();
                         strands=BothStrands())
threshold_sites = selectsites(pwm, sequences, ThresholdHits(5.0f0))

pfm = reconstruct_pfm(
    pwm,
    sequences,
    TopFractionHits(0.05);
    strands=BothStrands(),
    pseudocount=0.25f0,
)
```

Selectors are `BestPerSequence()`, `ThresholdHits(score)`, and
`TopFractionHits(fraction[, base_selector])`.

### Build and apply null distributions

```julia
models = [pwm, readmodel("examples/gata2.meme"),
          readmodel("examples/foxa2.meme")]
relations = parse_group_relations(
    "groups.tsv";
    known_names=Set(modelname(model) for model in models),
)

null_result = build_null(models, relations;
    sequences=sequences,
    background=background,
    metric=:co,
)

savenull("output/null_bundle", null_result.distribution)
null_distribution = loadnull("output/null_bundle")
annotated = annotate_results([result], null_distribution;
                             effective_number_of_targets=length(models))
```

### Add a custom model

Custom model families use multiple dispatch; there is no handler registry. A
model without additional context implements only `modelname`, `motif_length`,
and `scan_pair_kernel!`:

```julia
import Mimosa

struct ConsensusModel <: Mimosa.AbstractMotifModel
    label::String
    pattern::Vector{UInt8}
end

Mimosa.modelname(model::ConsensusModel) = model.label
Mimosa.motif_length(model::ConsensusModel) = length(model.pattern)

function Mimosa.scan_pair_kernel!(
    forward_scores::AbstractVector{Float32},
    reverse_scores::AbstractVector{Float32},
    model::ConsensusModel,
    sequence::AbstractVector{UInt8},
    n_positions::Int,
)
    pattern = model.pattern
    reverse_pattern = UInt8[b == Mimosa.N_CODE ? b : 0x03 - b
                            for b in reverse(pattern)]

    @inbounds for position in 1:n_positions
        fscore = 0.0f0
        rscore = 0.0f0
        for offset in eachindex(pattern)
            base = sequence[position + offset - 1]
            fscore += base == pattern[offset]
            rscore += base == reverse_pattern[offset]
        end
        forward_scores[position] = fscore
        reverse_scores[position] = rscore
    end
    return forward_scores, reverse_scores
end

model = ConsensusModel("ACGT consensus", Mimosa.encode_sequence("ACGT"))
Mimosa.validate_model(model; capability=:compare)
scores = Mimosa.scan(model, sequences; strands=Mimosa.BestStrand())
result = Mimosa.compare(model, pwm, sequences; metric=:co)
```

Models that read flanking bases additionally implement `left_context` and/or
`right_context`; both default to zero. Mimosa derives `window_size`, scan-track
length, and site offset from this geometry. Cache and null construction also
require a stable SHA-256 `model_fingerprint` method. Built-in portable bundles
currently store built-in model kinds only, so external packages should parse
their own formats and construct their model type directly.

See [Extending Mimosa](docs/src/extending.md) and the downstream contract test
for the complete extension and testing requirements.

## Parallel execution

Threading is opt-in twice: start Julia with multiple runtime threads and pass a
`ThreadedExecution` policy to the API.

```bash
JULIA_NUM_THREADS=4 julia --project=Mimosa.jl
```

```julia
scores = scan(pwm, sequences;
              strands=BestStrand(),
              execution=ThreadedExecution(4))
results = compare(query, targets;
                  metric=:co,
                  execution=ThreadedExecution(4))
```

`SerialExecution()` is always the API default. The CLI `--threads=N` selects
the same policy but cannot create Julia runtime threads and rejects values
larger than `Threads.nthreads()`. Parallel work occurs only at the highest
independent level and preserves result order and exact discrete fields.

## Storage and security

User-facing model and null storage uses bounded TOML manifests plus
checksum-verified binary blobs. Mimosa does not load pickle, joblib, Julia
`Serialization`, or evaluate input. Readers reject path traversal, symlink
escape, malformed metadata, unsupported shapes and dtypes, non-finite data,
oversized declarations, and checksum mismatches before accepting a bundle.
Writes are staged and committed atomically.

Current portable format versions are model 2, null 3, and cache 2. Null bundles
use the profile strategy only.

## Development

Run commands from the repository root:

```bash
julia --project=Mimosa.jl -e 'using Pkg; Pkg.test()'

julia --project=Mimosa.jl/test -e \
  'using JuliaFormatter; @assert format("Mimosa.jl/src"; overwrite=false); @assert format("Mimosa.jl/test"; overwrite=false)'

julia --project=Mimosa.jl/docs Mimosa.jl/docs/make.jl

JULIA_NUM_THREADS=4 julia --project=Mimosa.jl/benchmark \
  Mimosa.jl/benchmark/runbenchmarks.jl
```

See the [quick start](docs/src/quickstart.md), [CLI guide](docs/src/cli.md),
[API reference](docs/src/api.md), [model guide](docs/src/models.md),
[numerical compatibility contract](docs/src/numerical_compatibility.md), and
[architecture](docs/src/architecture.md).

## References

1. Gupta S. et al. (2007). Quantifying similarity between motifs. *Genome
   Biology*, 8, R24. <https://doi.org/10.1186/gb-2007-8-2-r24>
2. Lambert S. A. et al. (2016). Motif comparison based on similarity of binding
   affinity profiles. *Bioinformatics*, 32(22), 3504-3506.
   <https://doi.org/10.1093/bioinformatics/btw489>
3. Siebert M. and Soding J. (2016). Bayesian Markov models consistently
   outperform PWMs at predicting motifs in nucleotide sequences. *Nucleic
   Acids Research*, 44(13), 6055-6069.
   <https://doi.org/10.1093/nar/gkw521>

## License

MIT. See [LICENSE](LICENSE).
