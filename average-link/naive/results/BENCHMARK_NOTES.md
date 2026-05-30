# Naive parallel average-link — benchmark results

Quick summary of the first sweep on Salle info.

## Setup

- `hac_naive` vs `hac_seq` on synthetic 2-D clustered data (k=5, seed 42).
- N ∈ {200, 500, 1000, 1500, 2000, 3000}, threads ∈ {1, 2, 3, 4, 8}.
- 5 timed runs + 1 warmup per config, median reported.
- Hardware: Intel Core i7-13700K (16 cores / 24 threads), 124 GiB RAM, Linux 5.14.

## Speedup vs sequential (median)

|   N  | t=1  | t=2  | t=3  | t=4  | t=8  |
|-----:|:----:|:----:|:----:|:----:|:----:|
|  200 | 0.42 | 0.25 | 0.36 | 0.18 | 0.19 |
|  500 | 0.94 | 0.99 | 0.94 | 0.89 | 0.49 |
| 1000 | 0.86 | 1.21 | 1.42 | 1.52 | 1.21 |
| 1500 | 0.86 | 1.28 | 1.61 | 1.82 | 1.98 |
| 2000 | 0.88 | 1.33 | 1.70 | 1.96 | 2.70 |
| 3000 | 0.95 | 1.37 | 1.76 | 2.08 | 3.59 |

## Takeaways

- **Small N is overhead-bound.** Below N=500 the parallel version loses
  to sequential — thread spawn/join per merge dominates the work.
- **Crossover ~N=1000.** First clear wins for t≥3.
- **Best speedup observed: 3.59× at N=3000, t=8.** Sublinear because of
  the serial Phase B, the per-iteration reduction, and the spawn-per-merge
  cost.
- **No single best thread count** — t=4 wins at mid N, t=8 wins at large N.

Raw data in `naive_benchmark.csv`.
