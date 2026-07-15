# Release and Validation

## Platform Contract

Julia 1.10 is the minimum supported version. Development and current local
validation also use Julia 1.12. The package is pure Julia and has no Python
runtime dependency. Platform support is determined by the repository CI matrix;
do not infer support from historical release plans.

The package is not currently registered in General. Develop it from a local
clone with `Pkg.develop(path="/path/to/Mimosa.jl")`.

## Required Validation

From the repository root:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'

julia --project=test -e \
  'using JuliaFormatter; @assert format("src"; overwrite=false); @assert format("test"; overwrite=false)'

julia --project=docs docs/make.jl

julia --project=test/downstream test/downstream/runtests.jl
```

Threaded validation must start Julia with multiple runtime threads and pass an
explicit `ThreadedExecution` policy in the workflow being measured.

## Release Invariants

- `using Mimosa` performs no filesystem I/O, launches no work, prints nothing,
  and changes no global thread or BLAS setting.
- CLI JSON is stdout-only; diagnostics are stderr-only; exit codes are stable.
- Public storage uses model format 2, null format 3, cache format 2, and
  annotated-result schema 1.
- Bundle reads remain bounded and checksum verified; writes remain atomic.
- Numerical tolerances, security tests, Aqua/JET checks, and frozen fixtures are
  not weakened to obtain a green release.
- Direct matrix comparison and the `motif` command are not compatibility
  obligations.

## Benchmark Reporting

Warm compilation first, record Julia/CPU/runtime thread count, use identical
inputs and seed, and report median time plus allocations. Distinguish raw
profile comparison, prepared alignment, model scan-to-compare, and one-to-many
workloads. Threaded benchmarks must pass `ThreadedExecution(...)` explicitly.
