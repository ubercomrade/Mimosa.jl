# Reproducibility

## Determinism

### Within one Julia version

All Mimosa.jl operations are deterministic:
- Serial and threaded scanning produce **identical** results regardless of
  thread count (results written to pre-allocated slots indexed by position)
- Tie-breaking is deterministic: score, contributing site count, smaller
  absolute shift, then orientation priority `++ > +- > -+ > --`; complete
  ties retain the first visited shift
- `build_null` uses a fixed eligible-pair schedule independent of RNG
- No `push!` to shared arrays from multiple threads

### Across thread counts

`SerialExecution` and `ThreadedExecution(n)` produce **identical** results for
all model families and null distributions. This is verified by property tests
that compare serial vs threaded execution at 1, 2, and 4 threads.

### Across Julia versions

Float arithmetic is IEEE 754 and deterministic across Julia versions for the
same operations. However:
- `hash()` is **not stable** across Julia versions — cache keys use SHA-256
  content fingerprints instead
- Compiler optimizations may change instruction ordering for Float32
  operations, but results remain within documented tolerances

### Across languages

**Not bit-compatible** with Python/NumPy due to:
- Different accumulation order (NumPy pairwise summation vs Julia loop)
- Different RNG algorithm (MersenneTwister vs PCG64)
- Different GEV optimizer (native BFGS vs SciPy L-BFGS-B)

Results match within documented tolerance classes (see
[Numerical Compatibility](@ref)).

## RNG behavior

### Library functions

No library function uses the global RNG. `make_random_sequences` accepts an
explicit `seed` parameter and uses `MersenneTwister(seed)`.

### `build_null`

`build_null` does **not** generate random data — it uses the input models
directly for comparison. The eligible-pair schedule is deterministic based on
group structure. No `AbstractRNG` is needed for the basic null build.

### Future: random null sampling

If random sequence generation is added for null distribution sampling, the
API will accept `AbstractRNG` with stable seed derivation per task/pair,
independent of thread count.

## Cache stability

Cache keys are content-based SHA-256 fingerprints incorporating:
- Format schema version
- Algorithm name and version tag
- Model content (weights/representation bytes)
- Sequence content (encoded bytes)
- Configuration parameters

Keys are **stable across Julia sessions and versions**. Corrupted or
partial cache files result in cache misses, not errors.
