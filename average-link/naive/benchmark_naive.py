"""
Benchmark the naive thread-parallel average-link HAC against the sequential
matrix baseline.

For each dataset size N in --sizes:
  1) Generate (or reuse) a synthetic clustered CSV of N points.
  2) Sanity-check correctness once: confirm hac_naive (at max threads)
     produces a dendrogram identical to hac_seq's.
  3) Time hac_seq once (warmup) + --repeats times -> take the median.
  4) For each thread count in --threads:
       Time hac_naive --repeats times -> take the median.
       Compute speedup = T_seq_median / T_naive_median.
  5) Write a tidy results CSV and print a summary table.

The script measures algorithm wall-clock time as reported by the binaries
themselves ("HAC done in X ms"), which excludes file I/O. Each timing run
is invoked in a fresh process, so there's no warmup bleed-through between
runs (we still discard one untimed warmup run per binary to absorb OS-level
filesystem caching of the input CSV).

Usage (run from repo root):
  python3 average-link/naive/benchmark_naive.py
  python3 average-link/naive/benchmark_naive.py --sizes 500 1000 2000 \\
      --threads 1 2 3 4 --repeats 5 --out results/naive_bench.csv
"""

from __future__ import annotations

import argparse
import csv
import statistics
import subprocess
import sys
from pathlib import Path
from typing import Iterable

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
NAIVE_DIR = Path(__file__).resolve().parent
SEQ_BIN = REPO_ROOT / "average-link" / "hac_seq"
NAIVE_BIN = REPO_ROOT / "average-link" / "hac_naive"
DATA_DIR = NAIVE_DIR / "bench_data"
TMP_DIR = NAIVE_DIR / "bench_tmp"

# Tolerance for matching merge distances when validating correctness.
# Same value used in scripts/validate_average_link.py.
DIST_TOL = 1e-5


# ─── Dataset generation ──────────────────────────────────────────────────────

def generate_dataset(n: int, k: int, seed: int, path: Path) -> None:
    """Generate a clustered 2-D dataset matching scripts/generate_data.py."""
    rng = np.random.default_rng(seed)
    centres = rng.uniform(0, 10, size=(k, 2))
    sizes = rng.multinomial(n, [1 / k] * k)
    points = np.vstack([
        rng.normal(loc=centre, scale=0.4, size=(size, 2))
        for centre, size in zip(centres, sizes)
    ])
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(path, points, delimiter=",", fmt="%.6f")


def ensure_dataset(n: int, k: int, seed: int) -> Path:
    """Return path to a dataset of N points, generating it on demand."""
    path = DATA_DIR / f"bench_n{n}_k{k}_s{seed}.csv"
    if not path.exists():
        generate_dataset(n, k, seed, path)
    return path


# ─── Binary execution ────────────────────────────────────────────────────────

def parse_time_ms(stdout: str) -> float:
    """Extract the 'HAC done in X ms' value printed by both binaries."""
    for line in stdout.splitlines():
        if line.startswith("HAC done in"):
            # Format: "HAC done in 12.345 ms  (... merges)"
            return float(line.split()[3])
    raise RuntimeError(f"Could not find 'HAC done in' line in:\n{stdout}")


def run_seq(input_csv: Path, output_csv: Path) -> float:
    result = subprocess.run(
        [str(SEQ_BIN), str(input_csv), str(output_csv)],
        check=True, capture_output=True, text=True,
    )
    return parse_time_ms(result.stdout)


def run_naive(input_csv: Path, output_csv: Path, n_threads: int) -> float:
    result = subprocess.run(
        [str(NAIVE_BIN), str(input_csv), str(output_csv), str(n_threads)],
        check=True, capture_output=True, text=True,
    )
    return parse_time_ms(result.stdout)


# ─── Correctness check ──────────────────────────────────────────────────────

def load_dendrogram_sorted(path: Path) -> list[tuple[float, int]]:
    """Load (dist, new_size) pairs, sorted — order-independent merge fingerprint."""
    pairs: list[tuple[float, int]] = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            pairs.append((float(row["dist"]), int(row["new_size"])))
    return sorted(pairs)


def dendrograms_match(seq_csv: Path, naive_csv: Path) -> bool:
    """True iff the two dendrograms agree on (dist, new_size) up to DIST_TOL."""
    a = load_dendrogram_sorted(seq_csv)
    b = load_dendrogram_sorted(naive_csv)
    if len(a) != len(b):
        return False
    for (d1, s1), (d2, s2) in zip(a, b):
        if abs(d1 - d2) > DIST_TOL or s1 != s2:
            return False
    return True


# ─── Benchmark driver ────────────────────────────────────────────────────────

def median_time(fn, repeats: int) -> tuple[float, list[float]]:
    """Discard one warmup run, return (median, all timings) of `repeats` runs."""
    fn()  # warmup
    samples = [fn() for _ in range(repeats)]
    return statistics.median(samples), samples


def benchmark(
    sizes: Iterable[int],
    threads: Iterable[int],
    repeats: int,
    k: int,
    seed: int,
    out_csv: Path,
) -> list[dict]:
    TMP_DIR.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []

    print(f"{'N':>6} {'config':>12} {'median_ms':>11} {'min_ms':>10} "
          f"{'max_ms':>10} {'speedup':>8}")
    print("-" * 62)

    for n in sizes:
        data_csv = ensure_dataset(n, k, seed)
        seq_out = TMP_DIR / f"seq_n{n}.csv"
        naive_out = TMP_DIR / f"naive_n{n}.csv"

        # ── Correctness check: max threads vs sequential ────────────────────
        max_t = max(threads)
        run_seq(data_csv, seq_out)
        run_naive(data_csv, naive_out, max_t)
        if not dendrograms_match(seq_out, naive_out):
            print(f"!!! CORRECTNESS FAIL at N={n} (threads={max_t}). Aborting.",
                  file=sys.stderr)
            sys.exit(2)

        # ── Sequential baseline timing ───────────────────────────────────────
        seq_median, seq_samples = median_time(
            lambda: run_seq(data_csv, seq_out), repeats
        )
        rows.append({
            "N": n, "config": "seq", "n_threads": 0,
            "median_ms": seq_median,
            "min_ms": min(seq_samples), "max_ms": max(seq_samples),
            "samples_ms": ";".join(f"{s:.3f}" for s in seq_samples),
            "speedup": 1.0,
        })
        print(f"{n:>6} {'seq':>12} {seq_median:>11.3f} "
              f"{min(seq_samples):>10.3f} {max(seq_samples):>10.3f} "
              f"{1.0:>8.2f}")

        # ── Parallel timings ─────────────────────────────────────────────────
        for t in threads:
            naive_median, naive_samples = median_time(
                lambda: run_naive(data_csv, naive_out, t), repeats
            )
            speedup = seq_median / naive_median if naive_median > 0 else 0.0
            rows.append({
                "N": n, "config": f"naive_t{t}", "n_threads": t,
                "median_ms": naive_median,
                "min_ms": min(naive_samples), "max_ms": max(naive_samples),
                "samples_ms": ";".join(f"{s:.3f}" for s in naive_samples),
                "speedup": speedup,
            })
            print(f"{n:>6} {f'naive t={t}':>12} {naive_median:>11.3f} "
                  f"{min(naive_samples):>10.3f} {max(naive_samples):>10.3f} "
                  f"{speedup:>8.2f}")

        print()

    # ── Write results CSV ────────────────────────────────────────────────────
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["N", "config", "n_threads", "median_ms",
                        "min_ms", "max_ms", "samples_ms", "speedup"],
        )
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} rows to {out_csv}")
    return rows


# ─── CLI ─────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sizes", type=int, nargs="+",
                        default=[200, 500, 1000, 1500, 2000],
                        help="Dataset sizes N to benchmark")
    parser.add_argument("--threads", type=int, nargs="+",
                        default=[1, 2, 3, 4, 8],
                        help="Thread counts to test (focus on 1-4 per project brief)")
    parser.add_argument("--repeats", type=int, default=5,
                        help="Timed runs per configuration (median is reported)")
    parser.add_argument("--k", type=int, default=5,
                        help="Number of underlying clusters in synthetic data")
    parser.add_argument("--seed", type=int, default=42,
                        help="RNG seed for synthetic data generation")
    parser.add_argument("--out", type=Path,
                        default=NAIVE_DIR / "results" / "naive_benchmark.csv",
                        help="Output CSV path for the results")
    args = parser.parse_args()

    # Sanity check binaries exist before doing anything else.
    for bin_path in (SEQ_BIN, NAIVE_BIN):
        if not bin_path.exists():
            print(f"ERROR: binary not found at {bin_path}", file=sys.stderr)
            print("Build it first: cd average-link && make seq naive",
                  file=sys.stderr)
            sys.exit(1)

    benchmark(
        sizes=args.sizes,
        threads=args.threads,
        repeats=args.repeats,
        k=args.k,
        seed=args.seed,
        out_csv=args.out,
    )


if __name__ == "__main__":
    main()
