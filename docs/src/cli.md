# CLI

The CLI is a thin adapter over the public API. Successful results are JSON on
stdout; diagnostics and errors go to stderr.

## Installation

With Julia 1.12 or newer, install the executable through the experimental Pkg
Apps interface:

```julia
using Pkg
Pkg.Apps.add(url="https://github.com/ubercomrade/Mimosa.jl.git")
```

After registration in General, `Pkg.Apps.add("Mimosa")` is sufficient. Ensure
that `~/.julia/bin` is on `PATH`. Examples include `mimosa -- --help` and
`mimosa -- profile ...`; the separator separates Pkg App arguments.

If Mimosa is installed as a regular Julia package, use `Pkg.add` and invoke the
CLI module through Julia:

```julia
using Pkg
Pkg.add(url="https://github.com/ubercomrade/Mimosa.jl.git")
```

```bash
julia -m Mimosa profile examples/pif4.meme examples/gata2.meme \
  --model1-type pwm --model2-type pwm \
  --fasta examples/foreground.fa --metric co
```

Standalone archives attached to GitHub Releases include Julia and do not
require a Julia installation. Extract the entire archive and invoke
`Mimosa/bin/mimosa` (`Mimosa/bin/mimosa.exe` on Windows).

## Commands

### `profile`

Compare model-derived or precomputed score profiles:

```bash
JULIA_NUM_THREADS=4 julia -m Mimosa \
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

Pass `--cache-dir .mimosa-cache` to `profile` to persist prepared profiles.
A cache hit skips scanning, empirical normalization, and anchor collection when
the model, sequences, background, and threshold are unchanged.

### `build-null`

```bash
julia -m Mimosa build-null motifs/ \
  --model-type pwm --groups groups.tsv --output output/null_bundle \
  --fasta examples/foreground.fa --metric co
```

The relation file is TSV/CSV with configurable motif-name and group columns.
Only cross-group eligible pairs are compared. `--strict` and
`--min-null-targets` control insufficient-target handling. The output is always
a version-4 profile null bundle. `--jobs` is a deprecated alias for `--threads`.

### `cache clear`

```bash
julia -m Mimosa \
  cache clear --cache-dir .mimosa-cache
```

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
