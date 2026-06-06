from __future__ import annotations

import argparse
import re
import subprocess
import sys

import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data"
RESULTS_DIR = REPO_ROOT / "results"
PLOTS_DIR = RESULTS_DIR / "plots" / "analysis"
BENCHMARKS_DIR = RESULTS_DIR / "benchmarks"
GENERATE_SCRIPT = REPO_ROOT / "scripts" / "generate_data.py"

DATASET_SIZES = [500, 1000, 2000, 5000]
THREAD_COUNTS = [1, 2, 4, 8]
N_CLUSTERS = 10
N_RUNS = 3

# Real 2-D clustering datasets from milaan9/Clustering-Datasets (folder
# "02. Synthetic"), pre-fetched into data/ and converted to headerless x,y CSV.
# Ordered ascending by point count.
REAL_DATASETS = [
    ("flame",       DATA_DIR / "real_flame.csv"),
    ("R15",         DATA_DIR / "real_R15.csv"),
    ("aggregation", DATA_DIR / "real_aggregation.csv"),
    ("D31",         DATA_DIR / "real_D31.csv"),
]

FAMILIES = {
    "single": {
        "title": "Single-Link",
        "prefix": "single_link",
        "binaries": {
            "seq": REPO_ROOT / "single-link" / "single_link_baseline",
            "mst": REPO_ROOT / "single-link" / "single_link_mst",
        },
        "parallel": ["mst"],
        "labels": {"mst": "Boruvka MST parallel"},
        "colors": {"mst": "forestgreen"},
        "time_title": "Single-link wall-clock time vs dataset size",
    },
    "average": {
        "title": "Average-Link",
        "prefix": "average_link",
        "binaries": {
            "seq": REPO_ROOT / "average-link" / "hac_seq",
            "naive": REPO_ROOT / "average-link" / "hac_naive",
            "ppop": REPO_ROOT / "average-link" / "hac_ppop",
        },
        "parallel": ["naive", "ppop"],
        "labels": {"naive": "Naive parallel", "ppop": "pPOP parallel"},
        "colors": {"naive": "steelblue", "ppop": "darkorange"},
        "time_title": "Average-link wall-clock time vs dataset size",
    },
}

TIME_PATTERNS = [
    re.compile(r"finished in\s+([\d.]+)\s*ms"),
    re.compile(r"done in\s+([\d.]+)\s*ms"),
]


def parse_time_ms(stdout: str) -> float | None:
    for pattern in TIME_PATTERNS:
        match = pattern.search(stdout)
        if match:
            return float(match.group(1))
    return None


def build_cmd(family: str, algo: str, input_csv: Path, output_csv: Path,
              n_threads: int | None = None) -> list[str]:
    binary = FAMILIES[family]["binaries"][algo]
    cmd = [str(binary), str(input_csv), str(output_csv)]
    if algo != "seq":
        cmd.append(str(n_threads))
    return cmd


def generate_dataset(n_points: int, seed: int = 42) -> Path:
    path = DATA_DIR / f"bench_{n_points}.csv"
    if not path.exists():
        subprocess.run(
            [sys.executable, str(GENERATE_SCRIPT),
             "--n", str(n_points), "--k", str(N_CLUSTERS),
             "--seed", str(seed), "--out", str(path)],
            check=True,
        )
        print(f"  Generated {path.name}")
    return path


def run_trials(family: str, algo: str, input_csv: Path, output_csv: Path,
               n_threads: int | None = None) -> float | None:
    cmd = build_cmd(family, algo, input_csv, output_csv, n_threads)
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


def benchmark_family(family: str) -> pd.DataFrame:
    cfg = FAMILIES[family]
    results_csv = BENCHMARKS_DIR / f"{cfg['prefix']}_results.csv"
    speedup_csv = BENCHMARKS_DIR / f"{cfg['prefix']}_speedup.csv"
    tmp_out = BENCHMARKS_DIR / f"{cfg['prefix']}_bench_tmp.csv"
    results = []

    print("=" * 60)
    print(f"HAC {cfg['title']} Benchmark")
    print("=" * 60)

    for n in DATASET_SIZES:
        print(f"\n── N={n} ──")
        input_csv = generate_dataset(n)

        print("  seq ...", end=" ", flush=True)
        ms = run_trials(family, "seq", input_csv, tmp_out)
        if ms is not None:
            print(f"{ms:.1f} ms")
            results.append({"algo": "seq", "n_points": n, "n_threads": 1, "time_ms": ms})
        else:
            print("FAILED")

        for algo in cfg["parallel"]:
            for t in THREAD_COUNTS:
                print(f"  {algo} t={t} ...", end=" ", flush=True)
                ms = run_trials(family, algo, input_csv, tmp_out, n_threads=t)
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
        return pd.DataFrame()

    df = pd.DataFrame(results)
    df.to_csv(results_csv, index=False)
    print(f"\nRaw results → {results_csv}")

    speedup_rows = []
    for algo in cfg["parallel"]:
        for n in DATASET_SIZES:
            seq_rows = df[(df["algo"] == "seq") & (df["n_points"] == n)]
            if seq_rows.empty:
                continue
            seq_ms = seq_rows["time_ms"].values[0]
            for t in THREAD_COUNTS:
                row = df[(df["algo"] == algo) & (df["n_points"] == n) &
                         (df["n_threads"] == t)]
                if row.empty:
                    continue
                ms = row["time_ms"].values[0]
                speedup_rows.append({
                    "algo": algo, "n_points": n, "n_threads": t,
                    "time_ms": ms, "speedup": seq_ms / ms,
                })

    df_sp = pd.DataFrame(speedup_rows)
    df_sp.to_csv(speedup_csv, index=False)
    print(f"Speedup results → {speedup_csv}")

    plot_family(family, df, df_sp)
    print_summary(family, df_sp)
    return df


def plot_family(family: str, df: pd.DataFrame, df_sp: pd.DataFrame) -> None:
    if df_sp.empty:
        return

    cfg = FAMILIES[family]
    prefix = cfg["prefix"]
    labels = cfg["labels"]
    colors = cfg["colors"]
    max_t = df_sp["n_threads"].max()

    for algo in df_sp["algo"].unique():
        sub = df_sp[df_sp["algo"] == algo]
        fig, ax = plt.subplots(figsize=(7, 5))
        for n in sorted(sub["n_points"].unique()):
            grp = sub[sub["n_points"] == n].sort_values("n_threads")
            ax.plot(grp["n_threads"], grp["speedup"], marker="o", label=f"N={n}")
        ax.plot([1, max_t], [1, max_t], "k--", linewidth=0.8, label="Ideal")
        ax.set_xlabel("Number of threads")
        ax.set_ylabel("Speedup vs sequential")
        ax.set_title(f"Speedup — {labels.get(algo, algo)}")
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_xticks(sorted(df_sp["n_threads"].unique()))
        fig.savefig(PLOTS_DIR / f"{prefix}_speedup_{algo}.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 5))
    for algo in df_sp["algo"].unique():
        sub = df_sp[(df_sp["algo"] == algo) &
                    (df_sp["n_threads"] == max_t)].sort_values("n_points")
        ax.plot(sub["n_points"], sub["speedup"], marker="o",
                color=colors.get(algo), label=labels.get(algo, algo))
    ax.axhline(1, color="k", linestyle="--", linewidth=0.8, label="Sequential")
    ax.set_xlabel("Dataset size (N points)")
    ax.set_ylabel(f"Speedup vs sequential (t={max_t})")
    ax.set_title(f"Speedup vs dataset size at t={max_t}")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.savefig(PLOTS_DIR / f"{prefix}_speedup_vs_size.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 5))
    sub = df[df["algo"] == "seq"].sort_values("n_points")
    if not sub.empty:
        ax.plot(sub["n_points"], sub["time_ms"], marker="o",
                color="gray", label="Sequential")
    for algo in cfg["parallel"]:
        sub = df[(df["algo"] == algo) & (df["n_threads"] == max_t)].sort_values("n_points")
        if not sub.empty:
            ax.plot(sub["n_points"], sub["time_ms"], marker="s",
                    color=colors.get(algo),
                    label=f"{labels.get(algo, algo)} (t={max_t})")
    ax.set_xlabel("Dataset size (N points)")
    ax.set_ylabel("Time (ms)")
    ax.set_title(cfg["time_title"])
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.savefig(PLOTS_DIR / f"{prefix}_time_vs_size.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Plots → {PLOTS_DIR}")


def print_summary(family: str, df_sp: pd.DataFrame) -> None:
    if df_sp.empty:
        return
    labels = FAMILIES[family]["labels"]
    print("\n" + "=" * 60)
    print(f"SPEEDUP vs sequential — {FAMILIES[family]['title']}")
    print("=" * 60)
    for algo in df_sp["algo"].unique():
        print(f"\n{labels.get(algo, algo)}:")
        pivot = df_sp[df_sp["algo"] == algo].pivot(
            index="n_threads", columns="n_points", values="speedup"
        )
        print(pivot.to_string(float_format=lambda x: f"{x:.2f}x"))


def count_points(csv_path: Path) -> int:
    """Number of non-empty rows in a headerless CSV."""
    with open(csv_path) as f:
        return sum(1 for line in f if line.strip())


def benchmark_family_real(family: str) -> pd.DataFrame:
    """Same as benchmark_family but iterates over the real datasets list."""
    cfg = FAMILIES[family]
    results_csv = BENCHMARKS_DIR / f"{cfg['prefix']}_real_results.csv"
    speedup_csv = BENCHMARKS_DIR / f"{cfg['prefix']}_real_speedup.csv"
    tmp_out = BENCHMARKS_DIR / f"{cfg['prefix']}_real_bench_tmp.csv"
    results = []

    print("=" * 60)
    print(f"HAC {cfg['title']} Benchmark — real datasets")
    print("=" * 60)

    for label, csv_path in REAL_DATASETS:
        if not csv_path.exists():
            print(f"\n!! {label}: {csv_path.name} missing, skipping")
            continue
        n = count_points(csv_path)
        print(f"\n── {label} (N={n}) ──")

        print("  seq ...", end=" ", flush=True)
        ms = run_trials(family, "seq", csv_path, tmp_out)
        if ms is not None:
            print(f"{ms:.1f} ms")
            results.append({"algo": "seq", "dataset": label, "n_points": n,
                            "n_threads": 1, "time_ms": ms})
        else:
            print("FAILED")

        for algo in cfg["parallel"]:
            for t in THREAD_COUNTS:
                print(f"  {algo} t={t} ...", end=" ", flush=True)
                ms = run_trials(family, algo, csv_path, tmp_out, n_threads=t)
                if ms is not None:
                    print(f"{ms:.1f} ms")
                    results.append({"algo": algo, "dataset": label,
                                    "n_points": n, "n_threads": t,
                                    "time_ms": ms})
                else:
                    print("FAILED")

    if tmp_out.exists():
        tmp_out.unlink()

    if not results:
        print("\nNo real-dataset results collected.")
        return pd.DataFrame()

    df = pd.DataFrame(results)
    df.to_csv(results_csv, index=False)
    print(f"\nRaw results → {results_csv}")

    speedup_rows = []
    for algo in cfg["parallel"]:
        for label, _ in REAL_DATASETS:
            seq_rows = df[(df["algo"] == "seq") & (df["dataset"] == label)]
            if seq_rows.empty:
                continue
            seq_ms = seq_rows["time_ms"].values[0]
            n_pts = seq_rows["n_points"].values[0]
            for t in THREAD_COUNTS:
                row = df[(df["algo"] == algo) & (df["dataset"] == label) &
                         (df["n_threads"] == t)]
                if row.empty:
                    continue
                ms = row["time_ms"].values[0]
                speedup_rows.append({
                    "algo": algo, "dataset": label, "n_points": n_pts,
                    "n_threads": t, "time_ms": ms, "speedup": seq_ms / ms,
                })

    df_sp = pd.DataFrame(speedup_rows)
    df_sp.to_csv(speedup_csv, index=False)
    print(f"Speedup results → {speedup_csv}")

    plot_family_real(family, df, df_sp)
    print_summary_real(family, df_sp)
    return df


def plot_family_real(family: str, df: pd.DataFrame, df_sp: pd.DataFrame) -> None:
    if df_sp.empty:
        return

    cfg = FAMILIES[family]
    prefix = cfg["prefix"]
    labels = cfg["labels"]
    colors = cfg["colors"]
    max_t = df_sp["n_threads"].max()

    dataset_order = (df_sp[["dataset", "n_points"]]
                     .drop_duplicates()
                     .sort_values("n_points")["dataset"].tolist())

    for algo in df_sp["algo"].unique():
        sub = df_sp[df_sp["algo"] == algo]
        fig, ax = plt.subplots(figsize=(7, 5))
        for label in dataset_order:
            grp = sub[sub["dataset"] == label].sort_values("n_threads")
            if grp.empty:
                continue
            n_pts = grp["n_points"].iloc[0]
            ax.plot(grp["n_threads"], grp["speedup"], marker="o",
                    label=f"{label} (N={n_pts})")
        ax.plot([1, max_t], [1, max_t], "k--", linewidth=0.8, label="Ideal")
        ax.set_xlabel("Number of threads")
        ax.set_ylabel("Speedup vs sequential")
        ax.set_title(f"Real datasets — {labels.get(algo, algo)}")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        ax.set_xticks(sorted(df_sp["n_threads"].unique()))
        fig.savefig(PLOTS_DIR / f"{prefix}_real_speedup_{algo}.png",
                    dpi=150, bbox_inches="tight")
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    seq_sub = (df[df["algo"] == "seq"].set_index("dataset")
               .reindex(dataset_order))
    if not seq_sub.empty:
        ax.plot(dataset_order, seq_sub["time_ms"], marker="o",
                color="gray", label="Sequential")
    for algo in cfg["parallel"]:
        sub = (df[(df["algo"] == algo) & (df["n_threads"] == max_t)]
               .set_index("dataset").reindex(dataset_order))
        if not sub.empty:
            ax.plot(dataset_order, sub["time_ms"], marker="s",
                    color=colors.get(algo),
                    label=f"{labels.get(algo, algo)} (t={max_t})")
    ax.set_yscale("log")
    ax.set_xlabel("Dataset")
    ax.set_ylabel("Time (ms, log scale)")
    ax.set_title(f"{cfg['title']} — real datasets, time per dataset")
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)
    fig.autofmt_xdate(rotation=30)
    fig.savefig(PLOTS_DIR / f"{prefix}_real_time_per_dataset.png",
                dpi=150, bbox_inches="tight")
    plt.close(fig)

    # Speedup vs dataset size at the largest thread count: one line per
    # parallel algorithm, x-axis is n_points so the size scaling story is
    # visible. Markers are annotated with the dataset name for context.
    fig, ax = plt.subplots(figsize=(8, 5))
    for algo in cfg["parallel"]:
        sub = (df_sp[(df_sp["algo"] == algo) &
                     (df_sp["n_threads"] == max_t)]
               .sort_values("n_points"))
        if sub.empty:
            continue
        ax.plot(sub["n_points"], sub["speedup"], marker="o",
                color=colors.get(algo),
                label=f"{labels.get(algo, algo)} (t={max_t})")
        # Annotate each point with the dataset name.
        for _, row in sub.iterrows():
            ax.annotate(row["dataset"],
                        (row["n_points"], row["speedup"]),
                        textcoords="offset points", xytext=(5, 5),
                        fontsize=8, color="0.3")
    ax.axhline(1, color="k", linestyle="--", linewidth=0.8, label="Sequential")
    ax.set_xlabel("Dataset size (N points)")
    ax.set_ylabel(f"Speedup vs sequential (t={max_t})")
    ax.set_title(f"{cfg['title']} — real datasets, speedup vs dataset size")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.savefig(PLOTS_DIR / f"{prefix}_real_speedup_vs_size.png",
                dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Real-dataset plots → {PLOTS_DIR}")


def print_summary_real(family: str, df_sp: pd.DataFrame) -> None:
    if df_sp.empty:
        return
    labels = FAMILIES[family]["labels"]
    print("\n" + "=" * 60)
    print(f"SPEEDUP — {FAMILIES[family]['title']} (real datasets)")
    print("=" * 60)
    for algo in df_sp["algo"].unique():
        print(f"\n{labels.get(algo, algo)}:")
        pivot = df_sp[df_sp["algo"] == algo].pivot(
            index="n_threads", columns="dataset", values="speedup"
        )
        print(pivot.to_string(float_format=lambda x: f"{x:.2f}x"))


def main() -> None:
    parser = argparse.ArgumentParser(description="Benchmark HAC implementations.")
    parser.add_argument(
        "--family", choices=["single", "average", "all"], default="all",
        help="Which linkage family to benchmark (default: all)",
    )
    parser.add_argument(
        "--datasets", choices=["synthetic", "real", "both"], default="synthetic",
        help="Which dataset suite to run (default: synthetic). 'real' uses "
             "the milaan9 datasets only, 'both' runs everything.",
    )
    args = parser.parse_args()

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    BENCHMARKS_DIR.mkdir(parents=True, exist_ok=True)

    families = ["single", "average"] if args.family == "all" else [args.family]
    for family in families:
        if args.datasets in ("synthetic", "both"):
            benchmark_family(family)
        if args.datasets in ("real", "both"):
            benchmark_family_real(family)


if __name__ == "__main__":
    main()
