"""
generate_data.py  –  create synthetic 2-D clustering datasets for benchmarking.

Usage:
    python scripts/generate_data.py --n 200 --k 5 --out data/test_200.csv
"""

import argparse
import numpy as np
from pathlib import Path


def generate(n: int, k: int, seed: int = 42) -> np.ndarray:
    rng = np.random.default_rng(seed)
    # Place k cluster centres randomly in [0, 10]^2
    centres = rng.uniform(0, 10, size=(k, 2))
    points = []
    sizes = rng.multinomial(n, [1 / k] * k)  # split n points across k clusters
    for c, s in zip(centres, sizes):
        pts = rng.normal(loc=c, scale=0.4, size=(s, 2))
        points.append(pts)
    return np.vstack(points)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n",   type=int, default=200,          help="Number of points")
    parser.add_argument("--k",   type=int, default=5,            help="Number of true clusters")
    parser.add_argument("--seed",type=int, default=42,           help="Random seed")
    parser.add_argument("--out", type=str, default="data/test.csv", help="Output CSV path")
    args = parser.parse_args()

    data = generate(args.n, args.k, args.seed)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(out, data, delimiter=",", fmt="%.6f")
    print(f"Saved {args.n} points ({args.k} clusters) → {out}")


if __name__ == "__main__":
    main()