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