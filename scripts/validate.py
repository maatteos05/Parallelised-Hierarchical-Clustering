from __future__ import annotations

import argparse
import csv
import math
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DENDRO_DIR = REPO_ROOT / "results" / "dendrograms"

TOL_SINGLE = 1e-6
TOL_AVERAGE = 1e-5

BINARIES = {
    "single_baseline": REPO_ROOT / "single-link" / "hac_single_baseline",
    "single_mst": REPO_ROOT / "single-link" / "hac_single_mst",
    "avg_seq": REPO_ROOT / "average-link" / "hac_seq",
    "avg_naive": REPO_ROOT / "average-link" / "hac_naive",
    "avg_ppop": REPO_ROOT / "average-link" / "hac_ppop",
}


def load_points(path: Path) -> list[list[float]]:
    points = []
    with path.open(newline="") as csv_file:
        for row in csv.reader(csv_file):
            if row:
                points.append([float(value) for value in row])
    return points


def load_dendrogram_pairs(path: Path) -> list[tuple[float, int]]:
    pairs = []
    with path.open(newline="") as csv_file:
        for row in csv.DictReader(csv_file):
            pairs.append((float(row["dist"]), int(row["new_size"])))
    return sorted(pairs)


def compare_pairs(
    reference: list[tuple[float, int]],
    candidate: list[tuple[float, int]],
    tol: float,
    ref_label: str,
    cand_label: str,
) -> bool:
    if len(reference) != len(candidate):
        print(f"  FAIL length mismatch: {ref_label}={len(reference)}, "
              f"{cand_label}={len(candidate)}")
        return False

    if all(abs(a - b) <= tol and sa == sb
           for (a, sa), (b, sb) in zip(reference, candidate)):
        print(f"  PASS {len(reference)} merges match (tol={tol})")
        return True

    for index, ((ref_dist, ref_size), (cand_dist, cand_size)) in enumerate(
            zip(reference, candidate)):
        if abs(ref_dist - cand_dist) > tol or ref_size != cand_size:
            print(f"  FAIL first mismatch at sorted row {index}:")
            print(f"    {ref_label}: dist={ref_dist:.8f}, size={ref_size}")
            print(f"    {cand_label}: dist={cand_dist:.8f}, size={cand_size}")
            break
    return False


def dist(a: list[float], b: list[float]) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def reference_average_link(points: list[list[float]]) -> list[tuple[float, int]]:
    n = len(points)
    max_c = 2 * n
    distances = [[0.0] * max_c for _ in range(max_c)]
    sizes = [0] * max_c

    for i in range(n):
        sizes[i] = 1
        for j in range(i + 1, n):
            d = dist(points[i], points[j])
            distances[i][j] = d
            distances[j][i] = d

    active = list(range(n))
    pairs = []
    next_id = n

    while len(active) > 1:
        best = math.inf
        c1 = c2 = -1
        for i, ci in enumerate(active):
            for cj in active[i + 1:]:
                if distances[ci][cj] < best:
                    best = distances[ci][cj]
                    c1, c2 = ci, cj

        new_id = next_id
        next_id += 1
        sizes[new_id] = sizes[c1] + sizes[c2]

        for other in active:
            if other in (c1, c2):
                continue
            w = (sizes[c1] * distances[c1][other]
                 + sizes[c2] * distances[c2][other]) / sizes[new_id]
            distances[new_id][other] = w
            distances[other][new_id] = w

        active = [c for c in active if c not in (c1, c2)] + [new_id]
        pairs.append((best, sizes[new_id]))

    return sorted(pairs)


def run_binary(binary: Path, input_csv: Path, output_csv: Path,
               threads: int | None = None) -> None:
    cmd = [str(binary), str(input_csv), str(output_csv)]
    if threads is not None:
        cmd.append(str(threads))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(result.stderr.strip() or f"{binary.name} failed")


def dendro_paths(stem: str) -> dict[str, Path]:
    return {
        "single_baseline": DENDRO_DIR / "single-link" / f"baseline_{stem}.csv",
        "single_mst": DENDRO_DIR / "single-link" / f"mst_{stem}.csv",
        "avg_seq": DENDRO_DIR / "average-link" / f"seq_{stem}.csv",
        "avg_naive": DENDRO_DIR / "average-link" / f"naive_{stem}.csv",
        "avg_ppop": DENDRO_DIR / "average-link" / f"ppop_{stem}.csv",
    }


def legacy_paths(stem: str) -> dict[str, Path]:
    """Fallback for dendrograms named e.g. baseline_100.csv for test_100.csv."""
    suffix = stem.removeprefix("test_") if stem.startswith("test_") else stem
    return {
        "single_baseline": DENDRO_DIR / "single-link" / f"baseline_{suffix}.csv",
        "single_mst": DENDRO_DIR / "single-link" / f"mst_{suffix}.csv",
        "avg_seq": DENDRO_DIR / "average-link" / f"seq_{suffix}.csv",
        "avg_naive": DENDRO_DIR / "average-link" / f"naive_{suffix}.csv",
        "avg_ppop": DENDRO_DIR / "average-link" / f"ppop_{suffix}.csv",
    }


def resolve_path(key: str, paths: dict[str, Path], legacy: dict[str, Path]) -> Path:
    path = paths[key]
    if path.exists():
        return path
    legacy_path = legacy[key]
    if legacy_path.exists():
        return legacy_path
    return path


def run_all_binaries(input_csv: Path, paths: dict[str, Path], threads: int) -> None:
    for key in paths:
        paths[key].parent.mkdir(parents=True, exist_ok=True)
    run_binary(BINARIES["single_baseline"], input_csv, paths["single_baseline"])
    run_binary(BINARIES["single_mst"], input_csv, paths["single_mst"], threads)
    run_binary(BINARIES["avg_seq"], input_csv, paths["avg_seq"])
    run_binary(BINARIES["avg_naive"], input_csv, paths["avg_naive"], threads)
    run_binary(BINARIES["avg_ppop"], input_csv, paths["avg_ppop"], threads)


def validate_single(paths: dict[str, Path], legacy: dict[str, Path]) -> bool:
    print("\n[single-link] baseline vs MST")
    baseline = resolve_path("single_baseline", paths, legacy)
    mst = resolve_path("single_mst", paths, legacy)
    if not baseline.exists() or not mst.exists():
        print(f"  SKIP missing files: {baseline.name}, {mst.name}")
        return False
    return compare_pairs(
        load_dendrogram_pairs(baseline),
        load_dendrogram_pairs(mst),
        TOL_SINGLE, "baseline", "mst",
    )


def validate_average(points: list[list[float]], paths: dict[str, Path],
                     legacy: dict[str, Path]) -> bool:
    reference = reference_average_link(points)
    ok = True
    for key, label in [("avg_seq", "seq"), ("avg_naive", "naive"), ("avg_ppop", "ppop")]:
        print(f"\n[average-link] {label} vs Python reference")
        path = resolve_path(key, paths, legacy)
        if not path.exists():
            print(f"  SKIP missing file: {path.name}")
            ok = False
            continue
        ok &= compare_pairs(
            reference,
            load_dendrogram_pairs(path),
            TOL_AVERAGE, "reference", label,
        )
    return ok


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate HAC implementations.")
    parser.add_argument("--input", required=True, help="Input point CSV")
    parser.add_argument("--threads", type=int, default=4,
                        help="Thread count for parallel binaries (default: 4)")
    parser.add_argument("--family", choices=["single", "average", "all"],
                        default="all")
    parser.add_argument("--run", action="store_true", default=True,
                        help="Run binaries before validating (default)")
    parser.add_argument("--no-run", action="store_false", dest="run",
                        help="Only validate existing dendrogram files")
    args = parser.parse_args()

    input_csv = Path(args.input)
    if not input_csv.is_file():
        print(f"Input not found: {input_csv}")
        sys.exit(1)

    stem = input_csv.stem
    paths = dendro_paths(stem)
    legacy = legacy_paths(stem)
    points = load_points(input_csv)

    print(f"Input: {input_csv} ({len(points)} points)")
    if args.run:
        print("Running binaries ...")
        try:
            run_all_binaries(input_csv, paths, args.threads)
        except RuntimeError as exc:
            print(f"ERROR: {exc}")
            sys.exit(1)

    ok = True
    if args.family in ("single", "all"):
        ok &= validate_single(paths, legacy)
    if args.family in ("average", "all"):
        ok &= validate_average(points, paths, legacy)

    print("\n" + ("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED"))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
