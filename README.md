# Parallelised Hierarchical Clustering

Team: Matteo Sainton, Mateo Fatas, Emeric Payer

This repository implements and benchmarks parallel hierarchical agglomerative clustering (HAC) algorithms:

- Single-link, parallelized through a minimum spanning tree (MST) / Boruvka-style construction, following Olson (1995) in a CPU multithreaded setting.
- Average-link, parallelized with a thread-parallel distance matrix approach, also following Olson (1995).
- Average-link with pPOP, using the partially overlapping partitioning scheme from Dash, Petrutiu and Scheuermann (2007).

Sequential baselines are provided for both linkage criteria. The implementations are in C++ and use `std::thread`-based primitives. Benchmarks use synthetic 2-D clustered datasets and measure runtime and speedup as thread count changes.

## Contributions

- Emeric Payer: sequential baseline and parallel single-link MST implementation.
- Mateo Fatas: naive parallel average-link implementation.
- Matteo Sainton: pPOP implementation on top of average-link.
- Common: data generation, validation, visualization, and benchmarking scripts.

## Repository Layout

| Path | Purpose |
|---|---|
| `data/` | Input point-cloud CSVs, one `x,y` pair per line, no header. |
| `results/` | Generated dendrograms, benchmark CSVs, plots, and visual outputs. |
| `scripts/` | Dataset generation, validation, dendrogram visualization, cluster plotting, and benchmark drivers. |
| `single-link/` | Single-link baseline and parallel MST implementation. |
| `average-link/` | Sequential, naive parallel, and pPOP average-link implementations. |

Data file patterns:

| Pattern | Description |
|---|---|
| `test*.csv` | Small correctness datasets, for example `test_100.csv`. |
| `bench_{N}.csv` | Benchmark inputs, with N = 500, 1000, 2000, 5000. |
| `bench_n*_k*_s*.csv` | Extra benchmark inputs from earlier local runs. |
| `sl_{N}_{run}.csv` | Single-link multi-seed inputs. |

Results layout:

| Folder | Contents |
|---|---|
| `results/dendrograms/single-link/` | Single-link merge CSVs in `cl1,cl2,dist,new_size` format. |
| `results/dendrograms/average-link/` | Average-link merge CSVs. |
| `results/benchmarks/` | Timing tables from `scripts/benchmark.py`. |
| `results/plots/analysis/` | Runtime and speedup plots from `scripts/benchmark.py`. |
| `results/plots/visualizations/` | Cluster colourings from `scripts/plot_clusters.py`. |
| `results/visual/` | ASCII dendrogram text from `scripts/visualize_dendrogram.py`. |

Legacy average-link benchmark material is also stored in `average-link/naive/results/`, especially `naive_benchmark.csv`.

## Build

From the repository root:

```sh
make
```

This builds:

- `average-link/hac_seq`: sequential average-link baseline.
- `average-link/hac_naive`: naive parallel average-link.
- `average-link/hac_ppop`: pPOP parallel average-link.
- `single-link/hac_single_baseline`: sequential single-link baseline.
- `single-link/hac_single_mst`: parallel single-link MST / Boruvka version.

Partial builds:

```sh
make average
make single
make clean
```

Inside `single-link/`, `make` also builds:

- `hac_single_baseline`
- `hac_single_mst`

## Quick Use

All commands in this section are run from the repository root.

Generate synthetic clustered 2-D data:

```sh
python scripts/generate_data.py --n 100 --k 4 --out data/test_100.csv
python scripts/generate_data.py --n 500 --k 10 --seed 42 --out data/bench_500.csv
```

Run single-link:

```sh
./single-link/hac_single_baseline \
  data/test_100.csv results/dendrograms/single-link/baseline_100.csv

./single-link/hac_single_mst \
  data/test_100.csv results/dendrograms/single-link/mst_100.csv 4
```

The last argument to `hac_single_mst` is the thread count.

Run average-link:

```sh
./average-link/hac_seq \
  data/test_100.csv results/dendrograms/average-link/seq_100.csv

./average-link/hac_naive \
  data/test_100.csv results/dendrograms/average-link/naive_100.csv 4

./average-link/hac_ppop \
  data/test_100.csv results/dendrograms/average-link/ppop_100.csv 4
```

Program arguments:

```sh
./single-link/hac_single_baseline <input.csv> <output.csv>
./single-link/hac_single_mst <input.csv> <output.csv> <threads>

./average-link/hac_seq <input.csv> <output.csv>
./average-link/hac_naive <input.csv> <output.csv> [threads]
./average-link/hac_ppop <input.csv> <output.csv> <threads> [n_cells_per_dim]
```

Each program prints wall-clock time to stdout.

If running from inside `single-link/` instead of the root:

```sh
./hac_single_baseline ../data/test_100.csv ../results/dendrograms/single-link/baseline_100.csv
./hac_single_mst ../data/test_100.csv ../results/dendrograms/single-link/mst_100.csv 4
python ../scripts/validate.py --no-run --input ../data/test_100.csv --family single
```

## Validation And Visualization

Check correctness against baselines and Python references:

```sh
python scripts/validate.py --input data/test_100.csv
```

Useful options:

- `--family single|average|all`
- `--threads 4`
- `--no-run`, to validate existing dendrograms without regenerating them

Render an ASCII dendrogram:

```sh
python scripts/visualize_dendrogram.py \
  results/dendrograms/average-link/ppop_100.csv average_ppop_100_vis.txt
```

The output is written to `results/visual/average_ppop_100_vis.txt`.

Plot cluster colourings by cutting a dendrogram after `n - k` merges:

```sh
python scripts/plot_clusters.py --input data/test_100.csv --k 4
python scripts/plot_clusters.py --input data/test_100.csv --k 4 --family average --run
```

Useful options:

- `--family single|average|all`
- `--threads 4`
- `--run`, to generate dendrograms first
- `--out <path>`
- `--dendro name:path`, repeatable

Cluster plots are written to `results/plots/visualizations/clusters_<stem>_k<k>_<family>.png`.

## Benchmarks

Run all benchmark families:

```sh
python scripts/benchmark.py
```

Run one family:

```sh
python scripts/benchmark.py --family single
python scripts/benchmark.py --family average
```

Benchmark scripts write CSVs to `results/benchmarks/` and figures to `results/plots/analysis/`. They require `pandas` and `matplotlib`.

Regenerate benchmark and cluster plots:

```sh
python scripts/benchmark.py
python scripts/plot_clusters.py --input data/test_100.csv --k 4
```

Main benchmark outputs:

- `results/benchmarks/single_link_results.csv`
- `results/benchmarks/average_link_results.csv`
- speedup tables and figures under `results/plots/analysis/`

## Implementation Notes

Single-link files:

- `single-link/baseline_olson_matrix.cpp`: sequential baseline following Olson's matrix-based idea instead of recomputing every point distance at every merge.
- `single-link/parallel_mst.cpp`: parallel version using the fact that single-link HAC can be recovered from an MST, built with a threaded Boruvka-style approach.

Single-link test flow:

```sh
python scripts/generate_data.py --n 100 --k 4 --out data/test_100.csv
./single-link/hac_single_baseline data/test_100.csv results/dendrograms/single-link/baseline_100.csv
./single-link/hac_single_mst data/test_100.csv results/dendrograms/single-link/mst_100.csv 4
python scripts/validate.py --no-run --input data/test_100.csv --family single
```



# TEMPORARY (CONTENT FROM THE OTHER MARKDOWN FILES!!!) TO REMOVE

## Recorded Single-Link Results

Correctness: all MST outputs matched the sequential baseline.

Speedup formula:

```text
speedup = baseline_time / mst_4_threads_time
```

Runtime results:

| n | seed | baseline (ms) | MST 4 threads (ms) | speedup |
|---:|---:|---:|---:|---:|
| 100 | 1 | 0.241 | 0.648 | 0.37x |
| 100 | 2 | 0.369 | 0.972 | 0.38x |
| 100 | 3 | 0.349 | 0.734 | 0.48x |
| 100 | 4 | 0.362 | 0.720 | 0.50x |
| 100 | 5 | 0.236 | 0.805 | 0.29x |
| 500 | 1 | 7.827 | 5.660 | 1.38x |
| 500 | 2 | 6.505 | 5.233 | 1.24x |
| 500 | 3 | 6.927 | 5.448 | 1.27x |
| 500 | 4 | 6.439 | 4.740 | 1.36x |
| 500 | 5 | 8.723 | 5.113 | 1.71x |
| 1000 | 1 | 29.153 | 19.043 | 1.53x |
| 1000 | 2 | 37.156 | 23.281 | 1.60x |
| 1000 | 3 | 38.086 | 28.730 | 1.33x |
| 1000 | 4 | 39.612 | 24.355 | 1.63x |
| 1000 | 5 | 33.365 | 18.216 | 1.83x |
| 2000 | 1 | 168.397 | 104.615 | 1.61x |
| 2000 | 2 | 170.629 | 107.922 | 1.58x |
| 2000 | 3 | 179.437 | 104.659 | 1.71x |
| 2000 | 4 | 138.731 | 101.477 | 1.37x |
| 2000 | 5 | 128.989 | 87.455 | 1.47x |

Averages:

| n | avg baseline (ms) | avg MST 4 threads (ms) | avg speedup |
|---:|---:|---:|---:|
| 100 | 0.311 | 0.776 | 0.40x |
| 500 | 7.284 | 5.239 | 1.39x |
| 1000 | 35.474 | 22.725 | 1.56x |
| 2000 | 157.237 | 101.226 | 1.55x |

Summary:

- For `n=100`, the MST version is slower because thread overhead dominates.
- From `n=500` onward, the 4-thread MST version is consistently faster than the sequential baseline.
- Observed average speedup for larger datasets is about `1.4x` to `1.6x`.

## Recorded Naive Average-Link Results

First sweep on Salle info:

- Compared `hac_naive` against `hac_seq` on synthetic 2-D clustered data with `k=5`, seed `42`.
- Used N in `{200, 500, 1000, 1500, 2000, 3000}` and threads in `{1, 2, 3, 4, 8}`.
- Reported the median of 5 timed runs plus 1 warmup per configuration.
- Hardware: Intel Core i7-13700K, 16 cores / 24 threads, 124 GiB RAM, Linux 5.14.

Speedup vs sequential:

| N | t=1 | t=2 | t=3 | t=4 | t=8 |
|---:|---:|---:|---:|---:|---:|
| 200 | 0.42 | 0.25 | 0.36 | 0.18 | 0.19 |
| 500 | 0.94 | 0.99 | 0.94 | 0.89 | 0.49 |
| 1000 | 0.86 | 1.21 | 1.42 | 1.52 | 1.21 |
| 1500 | 0.86 | 1.28 | 1.61 | 1.82 | 1.98 |
| 2000 | 0.88 | 1.33 | 1.70 | 1.96 | 2.70 |
| 3000 | 0.95 | 1.37 | 1.76 | 2.08 | 3.59 |

Takeaways:

- Small N is overhead-bound. Below N=500 the parallel version loses to sequential because thread spawn/join per merge dominates.
- Crossover is around N=1000, with the first clear wins for at least 3 threads.
- Best observed speedup is `3.59x` at N=3000 with 8 threads.
- Speedup is sublinear because of the serial Phase B, per-iteration reduction, and spawn-per-merge cost.
- No single thread count is best: 4 threads wins at mid N, 8 threads wins at large N.
- Raw data is in `average-link/naive/results/naive_benchmark.csv`.

## Recorded pPOP Average-Link Analysis

Final benchmark on an 8-core MacBook, speedup vs sequential baseline:

```text
pPOP:  N=500  N=1000  N=2000  N=5000
 t=1   1.16x  1.93x   4.01x   7.65x
 t=2   1.08x  2.02x   4.78x   9.22x
 t=4   0.99x  2.03x   5.02x  10.24x
 t=8   0.81x  2.21x   6.48x  17.42x
```

Paper Table 2 correspondence:

| Paper step | Our implementation | Match |
|---|---|---|
| 1-2. Each processor creates priority queues for its chunk of cells `Ccell(P)`, in parallel. | `build_heaps_thread` dispatched on the pool; threads claim cells via atomic `next_cell`. | Yes, parallel and dynamic. |
| 3. Repeat while `overall_min_dist < delta`. | Inner `while (n_active > 1)` with break on `best_dist >= delta`. | Yes. |
| 4-5. Each processor finds closest pair for each cell in its chunk. | `find_min_thread`, atomic cell counter, dynamic scheduling. | Yes, parallel and dynamic. |
| 6-7. Designated processor finds overall closest pair and its cell. | Sequential reduction over `results[]`. | Yes. |
| 8. Merge `CL1`, `CL2`. | Create `k`, set `sz[k]`, weighted `rep_x/rep_y`. | Yes. |
| 9-10. Each processor updates priority queues of `Cclus(C)`. | `update_heaps_thread` on the pool over `affected_clusters`. | Yes, parallel. |
| 11-12. Designated processor determines affected neighbouring cells. | Folded into `affected_clusters` build by scanning `k`'s cells. | Yes, sequential and cheap. |
| 13-15. Each processor updates priority queues of affected cells. | Same `update_heaps_thread` dispatch. | Yes, parallel. |

Partitioning follows POP: a `c x c` grid over the 2-D bounding box, with each cell expanded by `delta / 2` on every side so adjacent cells overlap by exactly `delta`. This guarantees that any two clusters closer than `delta` share at least one cell, so the merge is not missed. Empirically, every dendrogram matched `scipy` `linkage(method='average')` to `1e-5` across N=50..800 and t=1,4,8.

Nested control follows Section 2: `delta` starts at 0 and grows by x1.1 per outer iteration (`D_INCR=0.1`); `c` shrinks per iteration (`C_DECR=0.1`), the 90-10 rule. The final iteration, with small `c` and large `delta`, acts as the paper's Phase 2.

Important deviation: the paper's pPOP is built on centroid linkage, while this project uses average-link (UPGMA). For centroid linkage, a merged cluster has one representative point, so only the container cell must be updated. For average-link, there is no representative point. The merged distance must use the Lance-Williams combination:

```text
D[k][m] = (sz_i * D[i][m] + sz_j * D[j][m]) / (sz_i + sz_j)
```

Because this references distances to both parents, `D[k][m]` must be maintained globally for every active cluster `m`. This changes update work from `O(n/c)` to `O(n)` per merge before parallelization. We keep the paper's `O(N/(c*p))` advantage on find-min, but the update step is heavier: our global Lance-Williams update is `O(n/p)`, while the paper's container-only centroid update is `O((n/c) * log n / p)`. This is a linkage-choice cost, not a parallelization bug.

Per-merge complexity with `n` remaining clusters, `c` cells, and `p` processors:

| Step | Paper pPOP, centroid | Ours, average-link | Same? |
|---|---|---|---|
| Create priority queues per outer iteration | `O(N^2/(c*p))` | `O(N^2/(c*p))` | Yes |
| Find closest pair | `O(n/p)` | `O(n/p)` | Yes |
| Merge | `O(1)` | `O(1)` | Yes |
| Update queues / distances | `O((n/c) * log n / p)` | `O((n/p) * log n)` global | No, factor `c` lost |

Implementation findings:

- Global parallel Lance-Williams update is necessary for average-link correctness and makes the update `O(n/p)` instead of `O(n)`.
- A persistent thread pool is needed because the paper assumes persistent processors with negligible synchronization. Creating threads per merge broke that assumption.
- O(1) heap clear by swap-to-empty removed a large teardown cost: at N=2000, this dropped runtime from about 7100 ms to about 2745 ms.
- A cell-count decay floor fixed load balance. Naive x0.9 integer decay collapsed `c` from 14 to 1 in three iterations, leaving one cell with about 95% of merges. Keeping total cells around `2p` during the heavy-merge phase restored the paper's `c/p` cells-per-processor assumption and changed N=5000 from flat scaling around 14000 ms to 2848 -> 1251 ms.

How to read the pPOP speedups:

- Super-linear speedup vs sequential at N=5000 is expected because partitioning reduces total work and parallelism divides the reduced work.
- Scaling improves with N because larger datasets keep more clusters per cell for longer, so load balance holds longer.
- Small N remains weak because parallel work per iteration is tiny relative to synchronization and the single-cell endgame.

Limitations to state in the report:

- Benchmarks only go to N=5000 on a laptop. The stored matrix is `O((2N)^2)` doubles, about 800 MB at N=5000, and grows quadratically. N=30000 would need roughly 28 GB.
- The baseline is a lean flat-matrix scan, not the paper's pTRAD priority-queue baseline. Report speedup vs sequential and vs naive separately.
- The `O(c)` update reduction is unavailable for average-link. A centroid implementation could update only the container cell, but would use a different linkage criterion.