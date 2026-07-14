# Feature Matrix

| Area | Current support | Public API / command |
|---|---|---|
| Models | PWM (from MEME/PFM), BaMM, SiteGA, Dimont, Slim | `readmodel`, model constructors |
| Profiles | Precomputed ragged scores | `ScoreProfile`, `read_scores` |
| Sequences | FASTA and generated DNA | `readsequences`, `make_random_sequences` |
| Scanning | Forward, reverse, best, both | `scan`, `scan!`, `StrandPolicy` |
| Comparison | Profile-only scalar and one-to-many | `compare`, `prepare_profile` |
| Metrics | CO, row-wise CO, Dice, row-wise Dice, cosine | `AbstractProfileMetric` implementations |
| Sites | Best, threshold, top fraction | `selectsites` |
| Reconstruction | Orientation-aware PFM | `reconstruct_pfm` |
| Statistics | Native GEV, p-values, BH FDR, E-values | `fit_gev`, `pvalue`, `adjusted_pvalues`, `evalue` |
| Nulls | Profile strategy, format v3 | `build_null`, `savenull`, `loadnull` |
| Parallelism | Bounded deterministic tasks | `SerialExecution`, `ThreadedExecution` |
| Model storage | TOML/raw-binary format v2 | `writemodel`, `readmodel` |
| Cache | Explicit format-v2 content cache | `Cache`, `cache_get`, `cache_set`, `clearcache` |
| CLI | Five workflows | `profile`, `build-null`, `cache clear`, `inspect-model`, `convert-model` |

Direct matrix/tensor comparison, PCC/Euclidean motif metrics, the `motif` CLI
command, and the `"motif"` null strategy are deliberately removed.

GPU/distributed execution, ZIP bundles, automatic scientific-workflow cache
integration, and empirical GEV fallback are not current public contracts.
