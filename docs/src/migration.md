# Historical Python Migration

The repository now has one active implementation: Mimosa.jl. The former Python
package has been removed from the source tree. This page records the remaining
migration boundary for legacy data; it is not a dual-support policy.

## Architectural Changes

| Retired Python design | Current Julia design |
|---|---|
| generic models and string registries | concrete immutable model types and multiple dispatch |
| padded/masked batches | flat `EncodedSequenceBatch` and `RaggedArray` storage |
| Numba/global thread control | explicit `ExecutionPolicy` and bounded top-level tasks |
| direct matrix and profile strategies | profile-only comparison |
| pickle/joblib persistence | bounded TOML + checksum-verified model blobs |
| SciPy GEV fitting | native Float64 BFGS fitting |

Direct matrix alignment, PCC/Euclidean metrics, the `motif` CLI command, and
the `"motif"` null strategy were deliberately removed.

## Legacy Serialized Data

Mimosa.jl never reads pickle/joblib. Historical conversion helpers remain under
`scripts/` and must be used only with explicitly trusted input because
Python deserialization can execute arbitrary code. Convert outside the Julia
trust boundary, then validate the resulting portable bundle with `readmodel` or
`loadnull`.

Null bundles from the removed strategy or versions older than 2 are not accepted
by the current reader and must not be relabeled without recomputing a compatible
profile null distribution.

## Historical Fixtures

Package-local parser/model inputs remain under `test/fixtures/`. The
former root Python oracle corpus has been removed from the current source tree;
it is not evidence that retired Python APIs remain supported. Any replacement
corpus requires documented generators, dependency versions, commands,
checksums, and scientific review.
