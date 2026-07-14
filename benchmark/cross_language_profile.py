"""Reproducible Julia/Python profile benchmark on shared FASTA and MEME inputs."""

from __future__ import annotations

import argparse
import gc
import json
import os
import platform
import re
import subprocess
import tempfile
import time
from pathlib import Path

import numpy as np

from mimosa.comparison import compare_one_to_many, create_comparator_config
from mimosa.comparison.common import _select_best_orientation
from mimosa.comparison.profile import (
    PROFILE_ORIENTATION_PAIRS,
    _collect_anchor_sites,
    _score_profile_candidates,
)
from mimosa.functions import (
    prepare_profile_bundle,
    scores_to_empirical_log_tail_bundle,
)
from mimosa.io import read_fasta
from mimosa.models import read_models
from mimosa.scanning import scan_model_strands

N_SEQUENCES = 10_000
SEQUENCE_LENGTH = 100
N_TARGETS = 50
MOTIF_WIDTH = 15
SEED = 12_345
_RESULT_RE = re.compile(r"RESULT language=julia threads=(\d+) median_s=([0-9.]+) min_s=([0-9.]+)")
_STAGE_RE = re.compile(r"STAGE language=julia name=([^ ]+) median_s=([0-9.]+) min_s=([0-9.]+)")


def _write_inputs(directory: Path) -> tuple[Path, Path]:
    rng = np.random.default_rng(SEED)
    fasta = directory / "sequences.fa"
    motifs = directory / "motifs.meme"
    alphabet = np.asarray(list("ACGT"))

    with fasta.open("w", encoding="ascii") as handle:
        encoded = rng.integers(0, 4, size=(N_SEQUENCES, SEQUENCE_LENGTH), dtype=np.int8)
        for index, row in enumerate(encoded):
            handle.write(f">seq_{index}\n{''.join(alphabet[row])}\n")

    with motifs.open("w", encoding="ascii") as handle:
        handle.write("MEME version 4\n\nALPHABET= ACGT\n\n")
        handle.write("Background letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n")
        for index in range(N_TARGETS + 1):
            matrix = rng.dirichlet(np.full(4, 0.7), size=MOTIF_WIDTH)
            handle.write(f"MOTIF motif_{index:02d}\n")
            handle.write(f"letter-probability matrix: alength= 4 w= {MOTIF_WIDTH} nsites= 100\n")
            for row in matrix:
                handle.write("\t".join(f"{value:.9f}" for value in row) + "\n")
            handle.write("\n")
    return fasta, motifs


def _python_benchmark(fasta: Path, motifs: Path, threads: int, reps: int) -> dict[str, object]:
    sequences = read_fasta(fasta)
    models = read_models(motifs, "pwm")
    query, targets = models[0], models[1:]
    config = create_comparator_config(
        metric="co",
        n_jobs=threads,
        search_range=10,
        window_radius=5,
        realign_window=3,
        min_logfpr=0.0,
    )

    def workload():
        return compare_one_to_many(query, targets, "profile", config, sequences=sequences)

    raw_query = scan_model_strands(query, sequences)
    raw_targets = [scan_model_strands(target, sequences) for target in targets]
    normalized_query = scores_to_empirical_log_tail_bundle(raw_query)
    normalized_targets = [scores_to_empirical_log_tail_bundle(raw) for raw in raw_targets]
    prepared_query = prepare_profile_bundle(normalized_query)
    prepared_targets = [prepare_profile_bundle(bundle) for bundle in normalized_targets]

    def measure(stage):
        stage()
        times = []
        for _ in range(reps):
            gc.collect()
            started = time.perf_counter()
            stage()
            times.append(time.perf_counter() - started)
        return {"median_s": float(np.median(times)), "min_s": min(times)}

    def collect_anchors(bundle):
        threshold = None
        n_rows = int(bundle["values"].shape[1])
        return [
            _collect_anchor_sites(bundle["values"][strand], bundle["lengths"], threshold)
            for strand in range(bundle["values"].shape[0])
        ], n_rows

    stages = {
        "query_scan": measure(lambda: scan_model_strands(query, sequences)),
        "query_normalization": measure(lambda: scores_to_empirical_log_tail_bundle(raw_query)),
        "target_scan": measure(lambda: [scan_model_strands(target, sequences) for target in targets]),
        "target_normalization": measure(lambda: [scores_to_empirical_log_tail_bundle(raw) for raw in raw_targets]),
        "anchor_collection": measure(
            lambda: [collect_anchors(bundle) for bundle in [normalized_query, *normalized_targets]]
        ),
        "alignment_1v1": measure(
            lambda: _select_best_orientation(
                _score_profile_candidates(prepared_query, prepared_targets[0], PROFILE_ORIENTATION_PAIRS, config)
            )
        ),
        "prepared_1v50": measure(
            lambda: [
                _select_best_orientation(
                    _score_profile_candidates(prepared_query, target, PROFILE_ORIENTATION_PAIRS, config)
                )
                for target in prepared_targets
            ]
        ),
    }

    if len(workload()) != N_TARGETS:
        raise RuntimeError("Python warm-up did not return all targets")
    times = []
    for _ in range(reps):
        gc.collect()
        started = time.perf_counter()
        workload()
        times.append(time.perf_counter() - started)
    return {
        "language": "python",
        "threads": threads,
        "median_s": float(np.median(times)),
        "min_s": min(times),
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "stages": stages,
    }


def _julia_benchmark(repo: Path, fasta: Path, motifs: Path, threads: int, reps: int) -> dict[str, object]:
    runner = repo / "Mimosa.jl" / "benchmark" / "cross_language_profile.jl"
    env = os.environ.copy()
    env["JULIA_NUM_THREADS"] = str(threads)
    env.setdefault("JULIA_DEPOT_PATH", f"/tmp/mimosa-julia-depot:{Path.home() / '.julia'}")
    completed = subprocess.run(
        [
            "julia",
            f"--project={repo / 'Mimosa.jl' / 'benchmark'}",
            str(runner),
            str(fasta),
            str(motifs),
            str(threads),
            str(reps),
        ],
        check=True,
        capture_output=True,
        text=True,
        env=env,
    )
    match = _RESULT_RE.search(completed.stdout)
    if match is None:
        raise RuntimeError(f"Could not parse Julia result:\n{completed.stdout}\n{completed.stderr}")
    stages = {
        stage.group(1): {"median_s": float(stage.group(2)), "min_s": float(stage.group(3))}
        for stage in _STAGE_RE.finditer(completed.stdout)
    }
    return {
        "language": "julia",
        "threads": int(match.group(1)),
        "median_s": float(match.group(2)),
        "min_s": float(match.group(3)),
        "stages": stages,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threads", type=int, nargs="+", default=[1, 4])
    parser.add_argument("--reps", type=int, default=3)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.reps < 1 or any(thread < 1 for thread in args.threads):
        parser.error("reps and thread counts must be positive")

    repo = Path(__file__).resolve().parents[2]
    with tempfile.TemporaryDirectory(prefix="mimosa-cross-language-") as tmp:
        fasta, motifs = _write_inputs(Path(tmp))
        results = []
        for threads in args.threads:
            print(f"benchmarking Python with {threads} thread(s)...", flush=True)
            results.append(_python_benchmark(fasta, motifs, threads, args.reps))
            print(f"benchmarking Julia with {threads} thread(s)...", flush=True)
            results.append(_julia_benchmark(repo, fasta, motifs, threads, args.reps))

    for threads in args.threads:
        python_result = next(r for r in results if r["language"] == "python" and r["threads"] == threads)
        julia_result = next(r for r in results if r["language"] == "julia" and r["threads"] == threads)
        julia_result["speedup_vs_python"] = python_result["median_s"] / julia_result["median_s"]
    report = {
        "workload": {
            "mode": "profile",
            "sequences": N_SEQUENCES,
            "sequence_length": SEQUENCE_LENGTH,
            "query_motifs": 1,
            "target_motifs": N_TARGETS,
            "motif_width": MOTIF_WIDTH,
            "metric": "co",
            "timing_scope": "compute only; inputs loaded and JIT warm-up completed",
            "repetitions": args.reps,
            "seed": SEED,
        },
        "results": results,
    }
    rendered = json.dumps(report, indent=2)
    print(rendered)
    if args.output is not None:
        args.output.write_text(rendered + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
