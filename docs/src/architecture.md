# Architecture

## Design principles

1. **Library-first**: The public API is the Julia module; the CLI is a thin adapter
2. **Concrete domain types**: `PWM`, `BaMM`, `SiteGA`, `Dimont`, and `Slim` are parametric immutable
   structs — no `GenericModel` with `Any` fields
3. **Multiple dispatch**: Metrics, strand policies, execution policies are small
   types; `scan`, `compare`, `scorebounds` dispatch on model types
4. **Composition over hierarchy**: Independent aspects (strand policy, metric,
   execution) are separate types, not deep inheritance
5. **No string dispatch in hot paths**: String identifiers exist only at I/O/CLI
   boundary, converted to concrete types immediately
6. **Type stability**: All hot kernels return concrete types, zero per-position
   allocations
7. **Serial + top-level parallelism**: Inner kernels are serial and composable;
   parallelism is at the top level (sequences, targets, pairs)

## Module structure

```
src/
├── Mimosa.jl           # Module entry point, exports
├── errors.jl           # Error hierarchy
├── precompile.jl       # PrecompileTools workload
├── serialization.jl    # JSON serialization (to_json, to_dict)
├── cli.jl              # Thin CLI adapter
├── models/             # Model types and constructors
│   ├── models.jl       # Model hierarchy and includes
│   ├── pwm.jl          # PWM type, PFM conversion, geometry traits
│   ├── bamm.jl         # BaMM type
│   ├── sitega.jl       # SiteGA type
│   ├── dimont.jl       # Dimont type
│   ├── slim.jl         # Slim type
├── sequences/          # Sequence representation
│   ├── encoding.jl     # EncodedSequenceBatch, reverse_complement
│   ├── fasta.jl        # FASTA reader
│   └── ragged.jl       # RaggedArray
├── scanning/           # Scanning interface, validation, and kernels
│   ├── scanning.jl      # Scanner module includes
│   ├── strands.jl       # Strand policies
│   └── n_order_scan.jl # Checked boundary and shared rolling k-mer kernels
├── comparison/         # Profile comparison
│   ├── comparison.jl    # Comparison module includes
│   ├── results.jl      # ComparisonResult, compare()
│   └── profile_comparison.jl  # Profile alignment
├── profiles/           # Profile inputs, normalization, anchors, and alignment
│   ├── profiles.jl      # Profiles module includes
│   ├── precomputed_profile.jl # ScoreProfile precomputed profile source
│   ├── normalization.jl # EmpiricalLogTail
│   ├── anchors.jl       # Anchor collection
│   ├── metrics.jl       # Profile metrics
│   └── alignment.jl    # Shift-based window alignment
├── sites/              # Site extraction and PFM reconstruction
│   └── sites.jl
├── statistics/         # Null distributions and statistics
│   ├── statistics.jl    # Statistics module includes
│   ├── gev.jl          # Native GEV fit
│   ├── pvalues.jl      # p-value, BH FDR, E-value
│   ├── relations.jl    # Group relations
│   ├── null_distribution.jl  # NullDistribution, build_null
│   └── null_storage.jl # savenull, loadnull
├── io/                 # File format parsers
│   ├── io.jl            # I/O module includes and public readers
│   ├── pfm_readers.jl   # MEME, PFM parsers
│   ├── pwm_reader.jl    # PFM-to-PWM reader adapters
│   ├── bamm_reader.jl  # BaMM .ihbcp parser
│   ├── sitega_reader.jl # SiteGA .mat parser
│   ├── xml_parser.jl   # Minimal XML parser
│   ├── dimont_reader.jl # Dimont XML parser
│   ├── slim_reader.jl  # Slim XML parser
│   ├── score_reader.jl # Score profile reader
│   ├── bundle_storage.jl # Bounded bundle parsing helpers
│   └── model_storage.jl # Portable model bundle format
├── parallel/          # Execution policies
│   └── parallel.jl
└── cache/              # Content-based cache
    └── cache.jl
```

## ADRs

Architectural decisions are documented in the project's `docs/adr/` directory:

- ADR 0001: Model type hierarchy
- ADR 0002: Sequence representation
- ADR 0003: Storage format
- ADR 0004: Parallelism and RNG
- ADR 0005: GEV fitting
- ADR 0006: Coordinate/offset/orientation conventions

## Scanning contract

`scan(model, sequence)` returns a `Vector{Float32}` for forward, reverse, or
best-strand policies and a `StrandPair` for `BothStrands()`. Batch scans return
the corresponding `RaggedArray{Float32}` representation. `scan!` fills a
caller-provided destination for one sequence; it deliberately rejects
`BothStrands()` because that policy has two outputs. The shared kernels retain
serial floating-point reduction order, while batch-level work can be scheduled
through an explicit `ThreadedExecution` policy.

## Precompilation

`PrecompileTools` workload exercises representative paths during package
precompilation (not at `using Mimosa` time):
- PWM construction and scanning
- Motif comparison (all metrics)
- Site extraction and PFM reconstruction
- GEV fitting
- JSON serialization
- Higher-order model scanning
- Cache fingerprint computation
