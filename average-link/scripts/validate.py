"""
validate.py  –  compare the C++ HAC output against scipy's reference.

Usage:
    python scripts/validate.py --data data/test_200.csv --cpp_out data/dendrogram.csv

The script:
  1. Loads the raw data.
  2. Runs scipy average-link HAC (ground truth).
  3. Loads the C++ dendrogram CSV.
  4. Compares merge distances row-by-row (order must match; distances must
     agree to within a tight tolerance).
  5. Reports PASS / FAIL with details.
"""

import argparse
import sys
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage

TOL = 1e-5   # absolute tolerance for distance comparison
             # (1e-6 is too tight: Lance-Williams weighted averages accumulate
             #  floating-point error over many merge steps)


def run_scipy(data: np.ndarray) -> np.ndarray:
    """Return scipy linkage matrix (N-1 rows × 4 cols: cl1, cl2, dist, size)."""
    return linkage(data, method="average", metric="euclidean")


def load_cpp(path: str) -> np.ndarray:
    df = pd.read_csv(path)
    return df[["cl1", "cl2", "dist", "new_size"]].to_numpy(dtype=float)


def compare(scipy_dg: np.ndarray, cpp_dg: np.ndarray) -> bool:
    """
    Compare the two dendrograms.

    scipy and our C++ implementation may number merged clusters differently
    (scipy uses contiguous ids, ours starts from N), but the DISTANCES and
    SIZES at each step must be identical (assuming no ties that would allow
    different valid merge orders).

    Strategy: compare the sorted list of (dist, new_size) pairs.
    This is robust to different but equally-valid tie-breaking.
    """
    if len(scipy_dg) != len(cpp_dg):
        print(f"FAIL  length mismatch: scipy={len(scipy_dg)}, cpp={len(cpp_dg)}")
        return False

    # Extract (dist, size) pairs and sort
    scipy_pairs = np.sort(scipy_dg[:, 2:4], axis=0)   # sort by dist first
    cpp_pairs   = np.sort(cpp_dg[:, 2:4],   axis=0)

    dist_ok = np.allclose(scipy_pairs[:, 0], cpp_pairs[:, 0], atol=TOL, rtol=0)
    size_ok = np.array_equal(scipy_pairs[:, 1].astype(int),
                             cpp_pairs[:, 1].astype(int))

    if dist_ok and size_ok:
        print(f"PASS  {len(scipy_dg)} merges match (tol={TOL})")
        return True

    # Report first mismatch
    for i, (sp, cp) in enumerate(zip(scipy_pairs, cpp_pairs)):
        if abs(sp[0] - cp[0]) > TOL or int(sp[1]) != int(cp[1]):
            print(f"FAIL  first mismatch at sorted row {i}:")
            print(f"      scipy : dist={sp[0]:.8f}  size={int(sp[1])}")
            print(f"      cpp   : dist={cp[0]:.8f}  size={int(cp[1])}")
            break
    return False


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data",    required=True, help="Input CSV (raw data points)")
    parser.add_argument("--cpp_out", required=True, help="C++ dendrogram CSV")
    args = parser.parse_args()

    data     = np.loadtxt(args.data, delimiter=",")
    scipy_dg = run_scipy(data)
    cpp_dg   = load_cpp(args.cpp_out)

    print(f"Points : {len(data)}")
    print(f"Merges : scipy={len(scipy_dg)}  cpp={len(cpp_dg)}")
    print()

    ok = compare(scipy_dg, cpp_dg)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()