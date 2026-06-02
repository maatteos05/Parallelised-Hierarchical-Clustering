from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
PLOTS_DIR = REPO_ROOT / "results" / "plots" / "visualizations"

sys.path.insert(0, str(SCRIPT_DIR))
from validate import (  # noqa: E402
    dendro_paths,
    legacy_paths,
    load_points,
    resolve_path,
    run_all_binaries,
)

ALGO_LABELS = {
    "avg_seq": "Average-link (seq)",
    "avg_naive": "Average-link (naive)",
    "avg_ppop": "Average-link (pPOP)",
    "single_baseline": "Single-link (baseline)",
    "single_mst": "Single-link (MST)",
}

FAMILY_ALGOS = {
    "average": ["avg_seq", "avg_naive", "avg_ppop"],
    "single": ["single_baseline", "single_mst"],
    "all": ["avg_seq", "avg_naive", "avg_ppop", "single_baseline", "single_mst"],
}


def labels_from_dendrogram(dendro_path: Path, n: int, k: int) -> np.ndarray:
    """Return cluster label 0..k-1 for each point after cutting at k clusters."""
    if not 1 <= k <= n:
        raise ValueError(f"k must be between 1 and n={n}, got {k}")

    max_id = 2 * n
    parent = list(range(max_id))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    merges_to_apply = n - k
    with dendro_path.open(newline="") as csv_file:
        reader = csv.DictReader(csv_file)
        for merge_index, row in enumerate(reader):
            if merge_index >= merges_to_apply:
                break
            cl1 = int(row["cl1"])
            cl2 = int(row["cl2"])
            union(cl1, cl2)
            # New cluster IDs are assigned sequentially from n in merge order.
            union(n + merge_index, cl1)

    roots = [find(point_id) for point_id in range(n)]
    unique_roots = sorted(set(roots))
    root_to_label = {root: label for label, root in enumerate(unique_roots)}
    return np.array([root_to_label[root] for root in roots], dtype=int)


def parse_dendro_arg(spec: str) -> tuple[str, Path]:
    if ":" not in spec:
        raise argparse.ArgumentTypeError(
            f"Expected name:path, got {spec!r}"
        )
    name, path_str = spec.split(":", 1)
    name = name.strip()
    if not name:
        raise argparse.ArgumentTypeError(f"Empty name in {spec!r}")
    return name, Path(path_str)


def plot_panel(ax, xs: np.ndarray, ys: np.ndarray, labels: np.ndarray | None,
               title: str, k: int | None = None) -> None:
    if labels is None:
        ax.scatter(xs, ys, c="0.65", s=28, edgecolors="0.35", linewidths=0.4)
    else:
        cmap = plt.get_cmap("tab10" if k <= 10 else "tab20")
        ax.scatter(xs, ys, c=labels, cmap=cmap, s=28,
                   vmin=0, vmax=max(k - 1, 0), edgecolors="0.2", linewidths=0.3)
    ax.set_title(title, fontsize=10)
    ax.set_aspect("equal", adjustable="datalim")
    ax.grid(True, alpha=0.25, linewidth=0.5)


def default_output_path(input_csv: Path, k: int, family: str) -> Path:
    return PLOTS_DIR / f"clusters_{input_csv.stem}_k{k}_{family}.png"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot raw points and HAC cluster colourings at k clusters.",
    )
    parser.add_argument("--input", required=True, help="Input point CSV (x,y)")
    parser.add_argument("--k", type=int, required=True,
                        help="Target number of clusters (dendrogram cut)")
    parser.add_argument("--family", choices=["single", "average", "all"],
                        default="all", help="Which dendrograms to plot (default: all)")
    parser.add_argument("--dendro", action="append", default=[],
                        metavar="NAME:PATH",
                        help="Extra dendrogram to plot (repeatable)")
    parser.add_argument("--threads", type=int, default=4,
                        help="Thread count when running parallel binaries (default: 4)")
    parser.add_argument("--run", action="store_true", default=False,
                        help="Run HAC binaries before plotting")
    parser.add_argument("--out", type=str, default=None,
                        help="Output PNG path (default: results/plots/visualizations/clusters_<stem>_k<k>_<family>.png)")
    args = parser.parse_args()

    input_csv = Path(args.input)
    if not input_csv.is_file():
        print(f"Input not found: {input_csv}")
        sys.exit(1)

    points = load_points(input_csv)
    n = len(points)
    if args.k < 1 or args.k > n:
        print(f"ERROR: k must be between 1 and {n}, got {args.k}")
        sys.exit(1)

    xs = np.array([p[0] for p in points])
    ys = np.array([p[1] for p in points])

    stem = input_csv.stem
    paths = dendro_paths(stem)
    legacy = legacy_paths(stem)

    if args.run:
        print("Running HAC binaries ...")
        try:
            run_all_binaries(input_csv, paths, args.threads)
        except RuntimeError as exc:
            print(f"ERROR: {exc}")
            sys.exit(1)

    panels: list[tuple[str, Path | None]] = [("Raw points", None)]

    for key in FAMILY_ALGOS[args.family]:
        dendro_path = resolve_path(key, paths, legacy)
        if not dendro_path.is_file():
            print(f"WARNING: skipping {key}, missing {dendro_path}")
            continue
        panels.append((ALGO_LABELS[key], dendro_path))

    for name, dendro_path in args.dendro:
        if not dendro_path.is_file():
            print(f"WARNING: skipping custom {name}, missing {dendro_path}")
            continue
        panels.append((name, dendro_path))

    if len(panels) == 1:
        print("ERROR: no dendrogram files found. Run with --run or pass --dendro.")
        sys.exit(1)

    ncols = len(panels)
    fig, axes = plt.subplots(1, ncols, figsize=(4.2 * ncols, 4.2), squeeze=False)
    axes = axes[0]

    for ax, (title, dendro_path) in zip(axes, panels):
        if dendro_path is None:
            plot_panel(ax, xs, ys, None, title)
        else:
            labels = labels_from_dendrogram(dendro_path, n, args.k)
            plot_panel(ax, xs, ys, labels, f"{title}\n(k={args.k})", args.k)

    fig.suptitle(f"{input_csv.name} — {n} points, k={args.k}", fontsize=12, y=1.02)
    fig.tight_layout()

    out_path = Path(args.out) if args.out else default_output_path(input_csv, args.k, args.family)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved cluster plot to {out_path}")


if __name__ == "__main__":
    main()
