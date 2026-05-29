import argparse
import csv
import sys

TOL = 1e-6


def load_pairs(path: str) -> list[tuple[float, int]]:
    pairs: list[tuple[float, int]] = []
    with open(path, newline="") as csv_file:
        reader = csv.DictReader(csv_file)
        for row in reader:
            pairs.append((float(row["dist"]), int(row["new_size"])))
    return sorted(pairs)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seq", required=True, help="Sequential dendrogram CSV")
    parser.add_argument("--par", required=True, help="Parallel dendrogram CSV")
    args = parser.parse_args()

    seq = load_pairs(args.seq)
    par = load_pairs(args.par)

    if len(seq) != len(par):
        print(f"FAIL length mismatch: seq={len(seq)}, par={len(par)}")
        sys.exit(1)

    if all(abs(seq_dist - par_dist) <= TOL and seq_size == par_size
           for (seq_dist, seq_size), (par_dist, par_size) in zip(seq, par)):
        print(f"PASS {len(seq)} merges match (tol={TOL})")
        sys.exit(0)

    for index, (seq_row, par_row) in enumerate(zip(seq, par)):
        if abs(seq_row[0] - par_row[0]) > TOL or seq_row[1] != par_row[1]:
            print(f"FAIL first mismatch at sorted row {index}:")
            print(f"  seq: dist={seq_row[0]:.8f}, size={seq_row[1]}")
            print(f"  par: dist={par_row[0]:.8f}, size={par_row[1]}")
            break
    sys.exit(1)


if __name__ == "__main__":
    main()
