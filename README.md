# Mimosa.jl

Mimosa.jl is a Julia 1.12+ package and command-line tool for DNA motif comparison. It supports PWM/PFM, BaMM, SiteGA, Dimont, Slim, and precomputed score profiles. Mimosa compares heterogeneous motif models through their behavior on the same DNA sequences rather than by forcing their internal representations into a common matrix.

## Background

Transcription factors (TFs) are key regulators of gene expression. They bind specific DNA sequences, termed transcription factor binding sites (TFBSs), located in gene regulatory regions. Because TFBSs are highly variable, they are typically described using motifs, which capture the nucleotide preferences of a TF across a range of binding sites. The position weight matrix (PWM) is the de facto standard motif model, assuming independence between positions and additive nucleotide contributions. However, numerous studies have demonstrated the presence of dependencies between positions within TFBSs. To account for such dependencies, alternative models have been proposed, including Markov-based approaches (BaMM, InMoDe, Dimont, Slim, MODER2), discriminative methods (SiteGA), and deep learning models (DeepBind, DeepGRN, BERT-TFBS). However, a mature ecosystem of databases and tools has been developed primarily for PWMs, including motif comparison methods that are critical for result interpretation. Comparison of alternative motif models usually requires their conversion into PWMs, which inevitably leads to loss of information about positional dependencies. At present, no universal tool exists for their direct comparison. To address this gap, we developed MIMOSA (Model-Independent Motif Similarity Assessment), a tool that evaluates motif similarity based on their functional behavior, i.e., the scores they assign to DNA sequences rather than their internal model parameters.

Directly aligning the parameters of different model families is generally not meaningful. Mimosa instead uses the following profile comparison pipeline:

1. Scan sequences by motifs to get score profile.
2. Convert raw scores to empirical `-log10(FPR)` values, optionally using a
   separate background sequence set for calibration.
3. Select one best anchor per sequence when `-log10(FPR) == 0`, or all anchors
   above the requested threshold.
4. Compare site-centered windows over shifts and the four strand orientations.
5. Return the highest-scoring alignment, its offset, orientation, and number of
   contributing sites.

See [Method](docs/src/method.md) for metric definitions, statistical
significance, and references.

## Supported inputs

| Family | CLI type | Input |
|---|---|---|
| Precomputed profiles | `scores` | FASTA-like numeric score rows |
| PWM/PFM | `pwm` | MEME, plain PFM, or a Mimosa model bundle |
| BaMM | `bamm` | `.ihbcp` or a Mimosa model bundle |
| SiteGA | `sitega` | `.mat` or a Mimosa model bundle |
| Dimont | `dimont` | Jstacs XML or a Mimosa model bundle |
| Slim | `slim` | Jstacs XML or a Mimosa model bundle |

See [Supported Models](docs/src/models.md) for format and scanning details.

## Installation

Mimosa.jl requires Julia 1.12 or newer. Until the package is registered in the
General registry, install the API directly from the repository:

```julia
using Pkg
Pkg.add(url="https://github.com/ubercomrade/Mimosa.jl.git")
```

## CLI

Install the Julia-backed executable with the experimental Pkg Apps interface:

```julia
using Pkg
Pkg.Apps.add(url="https://github.com/ubercomrade/Mimosa.jl.git")
```

Add `~/.julia/bin` to `PATH`, then run `mimosa profile ...`. Use `--` before
application arguments when the first application argument starts with a dash,
for example `mimosa -- --help`.

The CLI can also be run directly from a repository checkout:

```bash
julia --project=. app/mimosa.jl \
  profile examples/pif4.meme examples/gata2.meme \
  --model1-type pwm --model2-type pwm \
  --fasta examples/foreground.fa --metric co
```

Successful commands write JSON to `stdout`; diagnostics go to `stderr`. See the
[CLI guide](docs/src/cli.md) for all commands, options, threading, null
distributions, and exit codes.

Standalone archives for Linux, Windows, and macOS are attached to tagged GitHub
Releases. They include the Julia runtime and do not require Julia to be
installed.

## Julia API

```julia
using Mimosa

query = readmodel("examples/pif4.meme")
target = readmodel("examples/gata2.meme")
sequences, _ = readsequences("examples/foreground.fa")

result = compare(query, target, sequences; metric=:co)
println(to_json(result))
```

See the [Quick Start](docs/src/quickstart.md) for scanning, prepared profiles,
site extraction, PFM reconstruction, null distributions, and threading. The
[API reference](docs/src/api.md) documents every public interface. Example
models, sequences, profiles, and runnable scripts are in
[`examples/`](examples/).

## Documentation

- [Method and similarity metrics](docs/src/method.md)
- [Quick Start](docs/src/quickstart.md)
- [CLI guide](docs/src/cli.md)
- [Julia API reference](docs/src/api.md)
- [Supported models and formats](docs/src/models.md)
- [Extending Mimosa with custom models](docs/src/extending.md)

## License

MIT. See [LICENSE](LICENSE).
