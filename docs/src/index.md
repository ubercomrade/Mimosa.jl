# Mimosa.jl

Mimosa.jl is a Julia 1.10+ package for DNA motif scanning, profile-based
comparison, site extraction, PFM reconstruction, and statistical evaluation.

It supports PWM models read from MEME or PFM input, BaMM, SiteGA, Dimont,
Slim, and precomputed [`ScoreProfile`](@ref) inputs. Model comparison is profile-only: models are
scanned against an explicit [`EncodedSequenceBatch`](@ref), empirically
normalized, anchored, and aligned across strands.

Supported comparison metrics are `co`, `co_rowwise`, `dice`, `dice_rowwise`,
and `cosine`. Direct matrix alignment, PCC/Euclidean motif metrics, the `motif`
CLI command, and the `"motif"` null strategy are not part of the current API.

## Capabilities

- Flat validated sequence and ragged score storage
- Serial and bounded threaded scanning for all model families
- Raw and prepared scalar/one-to-many profile comparison
- Best, threshold, and top-fraction site extraction
- Orientation-aware PFM reconstruction
- Native Float64 GEV fitting, p-values, BH FDR, and E-values
- Profile null construction and result annotation
- Versioned TOML/raw-binary model bundles and TOML/NPY null bundles with atomic writes
- Explicit content-addressed cache
- Thin JSON CLI with stable exit codes

## Package Overview

```@docs
Mimosa
```

## Navigation

- [Method](method.md)
- [Quick Start](quickstart.md)
- [Julia API](api.md)
- [CLI](cli.md)
- [Supported Models](models.md)
- [Feature Matrix](feature_matrix.md)
- [Data Layout](data_layout.md)
- [Numerical Compatibility](numerical_compatibility.md)
- [Reproducibility](reproducibility.md)
- [Storage Format](storage.md)
- [Security](security.md)
- [Historical Python Migration](migration.md)
- [Extending Mimosa](extending.md)
- [MotifHORDE Contract](downstream_contract.md)
- [Architecture](architecture.md)
- [Release and Validation](release.md)
