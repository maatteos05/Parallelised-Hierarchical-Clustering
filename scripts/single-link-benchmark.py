"""
single-link-benchmark.py — Benchmark single-link HAC implementations across dataset sizes and thread counts.

Run from the repo ROOT:
    python scripts/single-link-benchmark.py

Outputs:
    data/single_link_benchmark_results.csv   — raw timing data
    data/single_link_benchmark_speedup.csv   — speedup relative to sequential baseline
    data/plots/                              — speedup and time plots per algorithm / dataset size
"""

import subprocess
import re
import sys
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# ─────────────────────────────────────────────────────────────────────────────
# CONFIG
# ─────────────────────────────────────────────────────────────────────────────

REPO_ROOT       = Path(__file__).resolve().parent.parent
DATA_DIR        = REPO_ROOT / "data"
PLOTS_DIR       = DATA_DIR  / "plots"
RESULTS_CSV     = DATA_DIR  / "single_link_benchmark_results.csv"
SPEEDUP_CSV     = DATA_DIR  / "single_link_benchmark_speedup.csv"
GENERATE_SCRIPT = REPO_ROOT / "scripts" / "generate_data.py"

# Binaries — paths relative to REPO_ROOT
BINARIES = {
    "seq": REPO_ROOT / "single-link" / "hac_single_baseline",
    "mst": REPO_ROOT / "single-link" / "hac_single_mst",
}

# Dataset sizes to benchmark
DATASET_SIZES = [500, 1000, 2000, 5000]

# Thread counts for parallel implementations
THREAD_COUNTS = [1, 2, 4, 8]

# Number of clusters in generated data
N_CLUSTERS = 10

# Number of runs per configuration — take the minimum time for stability
N_RUNS = 3

# ─────────────────────────────────────────────────────────────────────────────
# CLI builders
# ─────────────────────────────────────────────────────────────────────────────

def build_cmd(algo, input_csv, output_csv, n_threads=None):
    if algo == "seq":
        return [str(BINARIES["seq"]), str(input_csv), str(output_csv)]
    elif algo == "mst":
        return [str(BINARIES["mst"]), str(input_csv), str(output_csv),
                str(n_threads)]
    raise ValueError(f"Unknown algorithm: {algo}")

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def parse_time_ms(stdout):
    """Extract time from a line containing 'finished in X ms'."""
    match = re.search(r"finished in\s+([\d.]+)\s*ms", stdout)
    return float(match.group(1)) if match else None


def generate_dataset(n_points, seed=42):
    """Generate a 2-D dataset if it does not already exist."""
    path = DATA_DIR / f"bench_{n_points}.csv"
    if not path.exists():
        subprocess.run(
            [sys.executable, str(GENERATE_SCRIPT),
             "--n", str(n_points),
             "--k", str(N_CLUSTERS),
             "--seed", str(seed),
             "--out", str(path)],
            check=True
        )
        print(f"  Generated {path.name}")
    return path


def run_trials(algo, input_csv, output_csv, n_threads=None):
    """Run N_RUNS trials and return the minimum time in ms, or None on failure."""
    cmd = build_cmd(algo, input_csv, output_csv, n_threads)
    times = []
    for _ in range(N_RUNS):
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
            if result.returncode != 0:
                print(f"    ERROR: {result.stderr.strip()}")
                return None
            ms = parse_time_ms(result.stdout)
            if ms is not None:
                times.append(ms)
        except subprocess.TimeoutExpired:
            print(f"    TIMEOUT ({algo} t={n_threads})")
            return None
        except FileNotFoundError:
            print(f"    BINARY NOT FOUND: {cmd[0]}")
            return None
    return min(times) if times else None

# ─────────────────────────────────────────────────────────────────────────────
# Main benchmark loop
# ─────────────────────────────────────────────────────────────────────────────

def main():
    DATA_DIR.mkdir(exist_ok=True)
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    results = []  # list of dicts: {algo, n_points, n_threads, time_ms}
    tmp_out = DATA_DIR / "single_link_bench_tmp.csv"

    print("=" * 60)
    print("HAC Single-Link Benchmark")
    print("=" * 60)

    for n in DATASET_SIZES:
        print(f"\n── N={n} ──")
        input_csv = generate_dataset(n)

        # Sequential baseline
        print(f"  seq ...", end=" ", flush=True)
        ms = run_trials("seq", input_csv, tmp_out)
        if ms is not None:
            print(f"{ms:.1f} ms")
            results.append({"algo": "seq", "n_points": n,
                             "n_threads": 1, "time_ms": ms})
        else:
            print("FAILED")

        # Parallel implementations
        for algo in ["mst"]:
            for t in THREAD_COUNTS:
                print(f"  {algo} t={t} ...", end=" ", flush=True)
                ms = run_trials(algo, input_csv, tmp_out, n_threads=t)
                if ms is not None:
                    print(f"{ms:.1f} ms")
                    results.append({"algo": algo, "n_points": n,
                                    "n_threads": t, "time_ms": ms})
                else:
                    print("FAILED")

    if tmp_out.exists():
        tmp_out.unlink()

    if not results:
        print("\nNo results collected — are the binaries compiled?")
        return

    # ── Save raw results ──────────────────────────────────────────────────────
    df = pd.DataFrame(results)
    df.to_csv(RESULTS_CSV, index=False)
    print(f"\nRaw results → {RESULTS_CSV}")

    # ── Compute speedup vs sequential ─────────────────────────────────────────
    speedup_rows = []
    for algo in ["mst"]:
        for n in DATASET_SIZES:
            seq_row = df[(df["algo"] == "seq") & (df["n_points"] == n)]
            if seq_row.empty:
                continue
            seq_ms = seq_row["time_ms"].values[0]

            for t in THREAD_COUNTS:
                row = df[(df["algo"] == algo) &
                         (df["n_points"] == n) &
                         (df["n_threads"] == t)]
                if row.empty:
                    continue
                ms = row["time_ms"].values[0]
                speedup_rows.append({
                    "algo":      algo,
                    "n_points":  n,
                    "n_threads": t,
                    "time_ms":   ms,
                    "speedup":   seq_ms / ms,
                })

    df_sp = pd.DataFrame(speedup_rows)
    df_sp.to_csv(SPEEDUP_CSV, index=False)
    print(f"Speedup results → {SPEEDUP_CSV}")

    # ── Plots ─────────────────────────────────────────────────────────────────
    plot_speedup(df_sp)
    plot_time_vs_size(df)
    print(f"Plots → {PLOTS_DIR}")

    print_summary(df_sp)


# ─────────────────────────────────────────────────────────────────────────────
# Plots
# ─────────────────────────────────────────────────────────────────────────────

ALGO_LABELS = {
    "mst": "Boruvka MST parallel",
}

COLORS = {
    "mst": "forestgreen",
}


def plot_speedup(df_sp):
    """One plot per algorithm: speedup vs threads, one line per dataset size."""
    if df_sp.empty:
        return

    for algo in df_sp["algo"].unique():
        sub = df_sp[df_sp["algo"] == algo]
        fig, ax = plt.subplots(figsize=(7, 5))

        for n in sorted(sub["n_points"].unique()):
            grp = sub[sub["n_points"] == n].sort_values("n_threads")
            ax.plot(grp["n_threads"], grp["speedup"],
                    marker="o", label=f"N={n}")

        max_t = df_sp["n_threads"].max()
        ax.plot([1, max_t], [1, max_t], "k--", linewidth=0.8, label="Ideal")

        ax.set_xlabel("Number of threads")
        ax.set_ylabel("Speedup vs sequential")
        ax.set_title(f"Speedup — {ALGO_LABELS.get(algo, algo)}")
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_xticks(sorted(df_sp["n_threads"].unique()))

        fig.savefig(PLOTS_DIR / f"single_link_speedup_{algo}.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

    # Combined speedup plot at max thread count
    if not df_sp.empty:
        max_t = df_sp["n_threads"].max()
        fig, ax = plt.subplots(figsize=(7, 5))
        for algo in df_sp["algo"].unique():
            sub = df_sp[(df_sp["algo"] == algo) &
                        (df_sp["n_threads"] == max_t)].sort_values("n_points")
            ax.plot(sub["n_points"], sub["speedup"],
                    marker="o", color=COLORS.get(algo),
                    label=ALGO_LABELS.get(algo, algo))
        ax.axhline(1, color="k", linestyle="--", linewidth=0.8, label="Sequential")
        ax.set_xlabel("Dataset size (N points)")
        ax.set_ylabel(f"Speedup vs sequential (t={max_t})")
        ax.set_title(f"Speedup vs dataset size at t={max_t}")
        ax.legend()
        ax.grid(True, alpha=0.3)
        fig.savefig(PLOTS_DIR / "single_link_speedup_vs_size.png", dpi=150, bbox_inches="tight")
        plt.close(fig)


def plot_time_vs_size(df):
    """Wall-clock time vs dataset size for all algorithms."""
    fig, ax = plt.subplots(figsize=(7, 5))
    max_t = max(THREAD_COUNTS)

    # Sequential
    sub = df[df["algo"] == "seq"].sort_values("n_points")
    if not sub.empty:
        ax.plot(sub["n_points"], sub["time_ms"],
                marker="o", color="gray", label="Sequential")

    # Parallel at max threads
    for algo in ["mst"]:
        sub = df[(df["algo"] == algo) &
                 (df["n_threads"] == max_t)].sort_values("n_points")
        if not sub.empty:
            ax.plot(sub["n_points"], sub["time_ms"],
                    marker="s", color=COLORS.get(algo),
                    label=f"{ALGO_LABELS.get(algo, algo)} (t={max_t})")

    ax.set_xlabel("Dataset size (N points)")
    ax.set_ylabel("Time (ms)")
    ax.set_title("Boruvka MST Wall-clock time vs dataset size")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.savefig(PLOTS_DIR / "single_link_time_vs_size.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def print_summary(df_sp):
    """Print speedup table to stdout."""
    if df_sp.empty:
        return
    print("\n" + "=" * 60)
    print("SPEEDUP vs sequential baseline")
    print("=" * 60)
    for algo in df_sp["algo"].unique():
        print(f"\n{ALGO_LABELS.get(algo, algo)}:")
        pivot = df_sp[df_sp["algo"] == algo].pivot(
            index="n_threads", columns="n_points", values="speedup"
        )
        print(pivot.to_string(float_format=lambda x: f"{x:.2f}x"))


if __name__ == "__main__":
    main()
