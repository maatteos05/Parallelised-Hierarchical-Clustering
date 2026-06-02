import argparse
from pathlib import Path

import numpy as np


def generate(n: int, k: int, seed: int = 42) -> np.ndarray:
    rng = np.random.default_rng(seed)
    centres = rng.uniform(0, 10, size=(k, 2))
    sizes = rng.multinomial(n, [1 / k] * k)

    points = []
    for centre, size in zip(centres, sizes):
        points.append(rng.normal(loc=centre, scale=0.4, size=(size, 2)))

    return np.vstack(points)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=200, help="Number of points")
    parser.add_argument("--k", type=int, default=5, help="Number of clusters")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--out", type=str, default="data/test.csv",
                        help="Output CSV path")
    args = parser.parse_args()

    data = generate(args.n, args.k, args.seed)
    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(output, data, delimiter=",", fmt="%.6f")
    print(f"Saved {args.n} points ({args.k} clusters) to {output}")


if __name__ == "__main__":
    main()
