# Mimosa.jl Benchmark Suite

Reproducible benchmark suite for Mimosa.jl, covering all representative workloads

## Quick start

```bash
# Run the full benchmark suite (prints human-readable + JSON to stdout)
julia --project=Mimosa.jl/benchmark Mimosa.jl/benchmark/runbenchmarks.jl

# Save results to a file
julia --project=Mimosa.jl/benchmark Mimosa.jl/benchmark/runbenchmarks.jl --output results.json

# Print environment metadata only
julia --project=Mimosa.jl/benchmark Mimosa.jl/benchmark/runbenchmarks.jl --report

# Compare against stored baseline
julia --project=Mimosa.jl/benchmark Mimosa.jl/benchmark/runbenchmarks.jl --baseline baseline.json

# Run the 1-vs-50 profile comparison benchmark
julia --project=Mimosa.jl/benchmark Mimosa.jl/benchmark/bench_1v50.jl

# Compare normalization designs on 500 × 20,000 calibration sequences
JULIA_NUM_THREADS=4 julia --project=Mimosa.jl/benchmark Mimosa.jl/benchmark/bench_normalization_options.jl

# Compare exact-tail and histogram-only normalization at increasing bin counts
JULIA_NUM_THREADS=4 julia --project=Mimosa.jl/benchmark Mimosa.jl/benchmark/bench_histogram_only_1v50.jl
```

## Thread configuration

Set `JULIA_NUM_THREADS` before launching Julia to control thread count:

```bash
JULIA_NUM_THREADS=4 julia --project=Mimosa.jl/benchmark Mimosa.jl/benchmark/runbenchmarks.jl
JULIA_NUM_THREADS=1 julia --project=Mimosa.jl/benchmark Mimosa.jl/benchmark/runbenchmarks.jl
```

## Options

| Option | Description |
|--------|-------------|
| `--output <file>` | Write JSON results to file (default: stdout) |
| `--report` | Print environment metadata only, no benchmarks |
| `--baseline <file>` | Compare results against stored baseline |
| `--samples <n>` | BenchmarkTools samples per benchmark (default: 100) |
| `--seconds <f>` | BenchmarkTools seconds budget per benchmark (default: 5.0) |
| `--help`, `-h` | Show help |

## Workloads covered

### PWM scanning
- Single-sequence forward, best-strand, and reverse complement
- Widths: 8, 15, 30
- Sequence lengths: 100, 200, 1000

### Batch scanning (ragged heavy-tail)
- Small (100 × 200), medium (1000 × 200), large (10000 × 200)
- Variable-length ragged batches (short-heavy, long-heavy)
- Serial and threaded scaling at 1/2/4 threads

### One-to-many profile comparison
- Target counts: 10, 100, 1000 (full suite)
- Dedicated 1-vs-50 benchmark (`bench_1v50.jl`): 10 000 × 100 bp sequences
- CO metric with search range and window radius

### Higher-order scanning (BaMM)
- Orders 0–5 (synthetic models)
- Real BaMM models from `examples/` (foxa2, gata2, gata4, myog — all order 4)

### Site extraction and PFM reconstruction
- BestPerSequence (low site density)
- ThresholdHits (high site density)
- TopFractionHits (medium density)
- Batch sizes: 100, 1000

### GEV fitting and statistics
- Sample sizes: 100, 500, 2000
- BH FDR correction (1000 p-values)

### Null distribution building
- Profile strategy with an explicit encoded sequence batch and CO metric
- Serial and threaded execution

### Storage (bundle write/read)
- PWM model round-trip
- BaMM model round-trip

### Startup and import latency
- `using Mimosa` in fresh subprocess
- `Pkg.precompile()` time
- CLI subprocess startup (`mimosa --version`)

## Metrics collected

Each benchmark records:

| Metric | Description |
|--------|-------------|
| `median_ns` | Median time in nanoseconds |
| `min_ns` | Minimum time in nanoseconds |
| `mean_ns` | Mean time in nanoseconds |
| `variance_ns` | Variance of time samples |
| `allocations` | Number of allocations |
| `memory_bytes` | Total bytes allocated |
| `n_samples` | Number of samples taken |
| `n_evals` | Number of evaluations per sample |
| `warmup` | Whether warm-up was performed (always true) |
| `parameters` | Workload-specific parameters (width, order, n_seqs, etc.) |

## Environment metadata

The report includes:

- Git commit SHA
- Julia version and executable path
- Machine architecture, OS, kernel
- CPU model, speed, cores
- Thread count
- Total and free RAM
- Package versions (Mimosa and dependencies)
- Timestamp
- Warm-up policy
- Sample count and seconds budget

## 1-vs-50 Normalization and Execution Benchmark

The `bench_1v50.jl` script compares two end-to-end execution strategies for
**1 query model vs 50 target models** on random DNA sequences.

### Configuration

| Parameter | Value |
|-----------|-------|
| Sequence length | 100 bp |
| Number of sequences | 10 000 |
| Number of targets | 50 |
| PWM width | 15 |
| Hybrid bins | 65 536 |
| Minimum log FPR | 3.0 |
| Metric | CO (overlap coefficient) |
| Search range | 10 |
| Window radius | 5 |
| Repetitions per measurement | 5 |

The benchmark compares two normalization strategies under the same
kernel-only parallel execution model:

| Case | Normalization | Target order | Computational kernels |
|------|---------------|--------------|-----------------------|
| Hybrid serial | `HybridEmpiricalLogTail` | Serial | `SerialExecution()` |
| Hybrid inner-threaded | `HybridEmpiricalLogTail` | Serial | `ThreadedExecution()` |
| Exact inner-threaded | `EmpiricalLogTail` | Serial | `ThreadedExecution()` |

Targets are processed one at a time. Scanning, normalization, anchor
collection, and profile alignment use the available Julia threads inside each
target. This bounds per-target working memory and removes nested execution
policies.

### Running

```bash
# One Julia thread (ThreadedExecution is effectively serial)
JULIA_NUM_THREADS=1 julia --project=Mimosa.jl/benchmark Mimosa.jl/benchmark/bench_1v50.jl

# Meaningful strategy comparison with four Julia threads
JULIA_NUM_THREADS=4 julia --project=Mimosa.jl/benchmark Mimosa.jl/benchmark/bench_1v50.jl
```

Optional workload overrides include `MIMOSA_BENCH_N_SEQUENCES`,
`MIMOSA_BENCH_N_TARGETS`, `MIMOSA_BENCH_N_REPS`, and
`MIMOSA_BENCH_HYBRID_BINS`. `MIMOSA_BENCH_MIN_LOGFPR` controls the exact-tail
cutoff and defaults to `3.0`; using `0.0` makes the Hybrid exact tail cover the
full calibration range. The report includes query preparation, 1-vs-1,
end-to-end 1-vs-50 timing, throughput, and result deltas between Hybrid and
Exact normalization.

The experimental `bench_histogram_only_1v50.jl` additionally removes the
exact-tail collection and sort entirely. It sweeps bin counts from the
comma-separated `MIMOSA_BENCH_HISTOGRAM_BINS` value (default:
`65536,262144,1048576`) and reports score error plus changed offsets,
orientations, and site counts relative to Exact normalization. Calibration
uses a separate 500 × 20 000 bp background by default; override it with
`MIMOSA_BENCH_BG_SEQUENCES` and `MIMOSA_BENCH_BG_LENGTH`.

### Historical comparison timings (obsolete benchmark design)

These retained results predate explicit normalization and execution strategy
selection and are not comparable with the current benchmark.

| Mode | 1 thread | 4 threads |
|------|-----------|-----------|
| 1-vs-1 (profiles precomputed) | 372 ms | 344 ms |
| **1-vs-50 (profiles precomputed)** | **19 301 ms** | **17 584 ms** |
| 1-vs-50 end-to-end (scan + compare) | 25 117 ms | 22 792 ms |

#### Historical throughput and per-target cost

| Metric | 1 thread | 4 threads |
|--------|-----------|-----------|
| Per-target (pure comparison) | 386 ms/target | 352 ms/target |
| Per-target (end-to-end) | 502 ms/target | 456 ms/target |
| Throughput (pure comparison) | 3 comparisons/sec | 3 comparisons/sec |
| Throughput (end-to-end) | 2 comparisons/sec | 2 comparisons/sec |
| Speedup (50 vs 50×1-vs-1, pure) | 1.0× | 1.0× |

### Historical conclusions (obsolete benchmark design)

1. **1-vs-50 comparison takes ~19 seconds** (pure comparison with precomputed
   profiles) and **~25 seconds** end-to-end (including target scanning). The
   throughput is **2–3 comparisons per second**.

2. **The bottleneck is profile preparation, not scanning or comparison.**
   Fitting the `EmpiricalLogTail` normalization over ~1.72 million scores
   (10 000 sequences × 86 positions × 2 strands) costs ~275 ms per profile —
   roughly **20× slower** than the PWM scan itself (~14 ms). The actual
   `profile_compare` call is negligible relative to preparation.

3. **Scanning is very fast.** A PWM scan over 10 000 × 100 bp sequences
   completes in ~14 ms. Precomputing all 50 target profiles (scan + profile
   build) takes ~920 ms total (~18 ms/model).

4. The historical ~1.1× result did not exercise the current kernel-level
   `ThreadedExecution`. Current one-to-many measurements use a single
   `execution=ThreadedExecution(n)` policy and record `Threads.nthreads()`.

5. **No scaling advantage from batching.** The 1-vs-50 time is exactly 50× the
   1-vs-1 time (speedup = 1.0×), confirming that each target is processed
   independently with no shared work reuse beyond the query profile.

### Historical sample comparison output

```
pwm_w15_s12446: score=0.4932 offset= 9 orient=++ n_sites=14952
pwm_w15_s12447: score=0.4977 offset=-1 orient=++ n_sites=16619
pwm_w15_s12448: score=0.4922 offset=-5 orient=++ n_sites=15654
pwm_w15_s12449: score=0.5381 offset= 0 orient=++ n_sites=16738
pwm_w15_s12450: score=0.5126 offset=-2 orient=++ n_sites=16720
```

## Regression baseline (E2)

The `baseline.json` file stores per-benchmark median timings for regression
comparison. Per `PLAN_2.md` E2:

- The baseline should be stored for a **controlled machine** or compared using
  stable normalized metrics.
- Scheduled CI publishes a comparison report but **does NOT block PRs** on noisy
  microbenchmarks.
- RC gate blocks only on **confirmed regressions** of agreed representative
  workloads.
- Any optimization change must include profile evidence and a compatibility rerun.

### Updating the baseline

```bash
# On a controlled machine with stable configuration:
julia --project=Mimosa.jl/benchmark Mimosa.jl/benchmark/runbenchmarks.jl \
  --output Mimosa.jl/benchmark/baseline.json
```

### Comparing against baseline

```bash
julia --project=Mimosa.jl/benchmark Mimosa.jl/benchmark/runbenchmarks.jl \
  --baseline Mimosa.jl/benchmark/baseline.json
```

Results with a ratio > 1.25× are flagged as potential regressions.

## Legacy benchmark file

The older `benchmarks.jl` file (Stage 9) is kept for backwards compatibility
but `runbenchmarks.jl` supersedes it with full PLAN_2.md E1/E2 coverage
including machine-readable JSON output, environment metadata, and baseline
comparison.

The `bench_1v50.jl` script is a standalone targeted benchmark for the
1-vs-50 profile comparison scenario. It is not part of the full
`runbenchmarks.jl` suite but can be run independently as described in the
[1-vs-50 Profile Comparison Benchmark](#1-vs-50-profile-comparison-benchmark)
section above.
