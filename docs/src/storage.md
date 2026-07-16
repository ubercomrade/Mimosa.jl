# Storage Format

Mimosa.jl stores models and null distributions as directory bundles. It never
uses pickle, joblib, or Julia `Serialization` for user-facing input.

## Model Bundles

Model bundles use format version 2:

```text
model_bundle/
├── manifest.toml
└── data/
    └── weights.bin  # PWM; higher-order models use representation.bin
```

The manifest records kind, name, dtype, shape, row-major layout, model-specific
geometry, exact payload length, relative blob path, and SHA-256 checksum. Data
blobs are raw little-endian `Float32` values. Supported kinds are `pwm`, `bamm`,
`sitega`, `dimont`, and `slim`. Version 1 NPY model bundles and PFM bundles are
legacy and rejected.

Use `writemodel(path, model)` and `readmodel(path)`.

## Null Bundles

Null bundles use format version 4:

```text
null_bundle/
├── manifest.toml
└── data/
    └── raw_null_scores.npy
```

The strategy is always `"profile"`. The manifest stores the profile metric,
Float64 GEV metadata, comparison counts, skipped queries, compatibility
fingerprints, and the checksum/shape declaration for raw scores.

Use `savenull(path, distribution)` and `loadnull(path)`.

## Binary Layout

In-memory matrices follow Julia column-major conventions. Model blobs are
explicitly row-major little-endian `Float32`; null-distribution bundles retain
their separate NPY format and Float64 payloads.

## Security and Atomicity

Readers reject traversal, absolute or escaping paths, symlink escape,
unsupported dtype/layout, shape or payload mismatch, non-finite values, invalid
checksums, and oversized declarations before allocation.

Writers build a complete sibling staging directory and commit it by rename.
Failed writes do not replace valid targets; orphan stages are ignored.

Any schema change requires a version increment, compatibility/migration checks,
tests, and updated documentation. Null-bundle format version 4 intentionally
rejects all older bundles; rebuild them before reuse.
