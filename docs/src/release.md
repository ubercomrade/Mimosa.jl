# Release and Validation

## Platform Contract

Julia 1.12 is the minimum supported version. The package is pure Julia and has
no Python runtime dependency. Platform support is determined by the repository
CI and release matrices.

The package is not currently registered in General. Develop it from a local
clone with `Pkg.develop(path="/path/to/Mimosa.jl")`.

## Distribution Channels

- General registration provides the Julia API through `Pkg.add("Mimosa")`.
- Pkg Apps provides the Julia-backed CLI through `Pkg.Apps.add("Mimosa")` and
  installs its launcher under `~/.julia/bin`.
- Tagged GitHub Releases provide PackageCompiler bundles that include Julia.

`CI.yml` tests the package on Linux, Windows, and macOS and verifies a real Pkg
Apps development installation. `Release.yml` accepts tags matching `v*`, first
validates that the tag equals the version in `Project.toml`, then tests and
builds platform-specific archives. A tag such as `v0.2.0` therefore requires
`version = "0.2.0"` in `Project.toml`.

`TagBot.yml` connects General registration to binary releases. Configure a
write-enabled SSH deploy key and store its private key in the repository secret
`TAGBOT_SSH_KEY`. TagBot then pushes the registered version tag through SSH,
which allows the tag event to trigger `Release.yml`. For the initial commit that
adds workflow files, create and push the tag manually if GitHub refuses the
automated release.

The release matrix currently publishes Linux x86-64/AArch64, Windows x86-64,
and macOS x86-64/AArch64 bundles. Each compiled executable is smoke-tested
before packaging; `SHA256SUMS` is generated when the GitHub Release is created.

## Required Validation

From the repository root:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'

julia --project=test -e \
  'using JuliaFormatter; @assert format("src"; overwrite=false); @assert format("test"; overwrite=false)'

julia --project=docs docs/make.jl

julia --project=test/downstream test/downstream/runtests.jl

julia --project=build build/build_app.jl
julia --project=build build/smoke_app.jl dist/Mimosa
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
