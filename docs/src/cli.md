# CLI

The CLI is a thin adapter over the public API. Successful results are JSON on
stdout; diagnostics and errors go to stderr.

Run commands from the repository root:

```bash
julia --project=Mimosa.jl Mimosa.jl/app/mimosa.jl <command> [options]
```

## Commands

### `profile`

Compare model-derived or precomputed score profiles:

```bash
JULIA_NUM_THREADS=4 julia --project=Mimosa.jl Mimosa.jl/app/mimosa.jl \
  profile examples/pif4.meme examples/gata2.meme \
  --model1-type pwm --model2-type pwm \
  --fasta examples/foreground.fa --metric co --threads 4
```

Required inputs are two positional model paths and `--model1-type` /
`--model2-type`. Types are `scores`, `pwm`, `bamm`, `sitega`, `dimont`, and
`slim`.

Comparison options include `--metric`, `--search-range`, `--window-radius`,
`--realign-window`, and `--min-logfpr`. Metric names are `co`, `co_rowwise`,
`dice`, `dice_rowwise`, and `cosine`.

Use `--fasta` for explicit scientific input or `--num-sequences`,
`--seq-length`, and `--seed` for generated sequences. `--background` accepts a
separate normalization FASTA. `--pvalue` requires a compatible explicit
`--null-distribution` bundle.

### `build-null`

```bash
julia --project=Mimosa.jl Mimosa.jl/app/mimosa.jl build-null motifs/ \
  --model-type pwm --groups groups.tsv --output output/null_bundle \
  --fasta examples/foreground.fa --metric co
```

The relation file is TSV/CSV with configurable motif-name and group columns.
Only cross-group eligible pairs are compared. `--strict` and
`--min-null-targets` control insufficient-target handling. The output is always
a version-3 profile null bundle. `--jobs` is a deprecated alias for `--threads`.

### `cache clear`

```bash
julia --project=Mimosa.jl Mimosa.jl/app/mimosa.jl \
  cache clear --cache-dir .mimosa-cache
```

### `inspect-model`

```bash
julia --project=Mimosa.jl Mimosa.jl/app/mimosa.jl \
  inspect-model examples/foxa2.ihbcp --type bamm
```

### `convert-model`

```bash
julia --project=Mimosa.jl Mimosa.jl/app/mimosa.jl \
  convert-model examples/pif4.meme output/pif4_bundle --type pwm
```

This converts a supported scientific model file into a portable Mimosa model
bundle. It does not deserialize pickle/joblib data.

## Threading

`--threads=N` selects `ThreadedExecution(N)` but cannot create Julia runtime
threads. Start Julia with `--threads=N` or `JULIA_NUM_THREADS=N`. The CLI rejects
a request larger than `Threads.nthreads()`.

## Process Contract

| Exit code | Meaning |
|---|---|
| 0 | success |
| 1 | usage or argument error |
| 2 | runtime, input, or scientific error |

Global flags are `--help`/`-h`, `--version`/`-V`, `--quiet`, and `--verbose`.
There are no interactive prompts.
