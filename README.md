# Parallelised Hierarchical Clustering

Team: Matteo Sainton, Mateo Fatas, Emeric Payer

This repository implements and benchmarks parallel hierarchical agglomerative clustering (HAC) algorithms:

- Single-link:
  - following Olson (1995) in a CPU multithreaded setting
  - parallelized through a minimum spanning tree (MST) / Boruvka-style construction
- Average-link:
  - parallelized with a thread-parallel distance matrix approach, also following Olson (1995)
  - with pPOP, using the partially overlapping partitioning scheme from Dash, Petrutiu and Scheuermann (2007)

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

## Build

From the repository root:

```sh
make
```

This builds:

- `average-link/hac_seq`: sequential average-link baseline.
- `average-link/hac_naive`: naive parallel average-link.
- `average-link/hac_ppop`: pPOP parallel average-link.
- `single-link/single_link_baseline`: sequential single-link baseline.
- `single-link/single_link_mst`: parallel single-link MST / Boruvka version.

Partial builds:

```sh
make average
make single
make clean
```

## Usage

All commands below are run from the repository root.

Use `DATASET_NAME` as the dataset file stem without `.csv`, for example `test_100`.

#### Generate synthetic clustered 2-D data:

```sh
python scripts/generate_data.py --n N_POINTS --k N_CLUSTERS --out data/DATASET_NAME.csv
python scripts/generate_data.py --n N_POINTS --k N_CLUSTERS --seed N_SEED --out data/DATASET_NAME.csv
```

Example:

```sh
python scripts/generate_data.py --n 100 --k 4 --out data/test_100.csv
python scripts/generate_data.py --n 100 --k 4 --seed 42 --out data/test_100.csv
```

#### Run single-link:

```sh
./single-link/single_link_baseline \
  data/DATASET_NAME.csv results/dendrograms/single-link/baseline_DATASET_NAME.csv

./single-link/single_link_mst \
  data/DATASET_NAME.csv results/dendrograms/single-link/mst_DATASET_NAME.csv N_THREADS
```

Example:

```sh
./single-link/single_link_baseline \
  data/test_100.csv results/dendrograms/single-link/baseline_test_100.csv

./single-link/single_link_mst \
  data/test_100.csv results/dendrograms/single-link/mst_test_100.csv 4
```

The last argument to `single_link_mst` is the thread count.

#### Run average-link:

```sh
./average-link/hac_seq \
  data/DATASET_NAME.csv results/dendrograms/average-link/seq_DATASET_NAME.csv

./average-link/hac_naive \
  data/DATASET_NAME.csv results/dendrograms/average-link/naive_DATASET_NAME.csv N_THREADS

./average-link/hac_ppop \
  data/DATASET_NAME.csv results/dendrograms/average-link/ppop_DATASET_NAME.csv N_THREADS
```

Example:

```sh
./average-link/hac_seq \
  data/test_100.csv results/dendrograms/average-link/seq_test_100.csv

./average-link/hac_naive \
  data/test_100.csv results/dendrograms/average-link/naive_test_100.csv 4

./average-link/hac_ppop \
  data/test_100.csv results/dendrograms/average-link/ppop_test_100.csv 4
```

#### Program arguments:

```sh
./single-link/single_link_baseline <input.csv> <output.csv>
./single-link/single_link_mst <input.csv> <output.csv> <threads>

./average-link/hac_seq <input.csv> <output.csv>
./average-link/hac_naive <input.csv> <output.csv> [threads]
./average-link/hac_ppop <input.csv> <output.csv> <threads> [n_cells_per_dim]
```

Each program prints wall-clock time to stdout.

## Validation And Visualization

#### Check correctness against baselines and Python references:

```sh
python scripts/validate.py --input data/DATASET_NAME.csv
```

Example:

```sh
python scripts/validate.py --input data/test_100.csv
```

Useful options:

- `--family single|average|all`
- `--threads N_THREADS`
- `--no-run`, to validate existing dendrograms without regenerating them

#### Render an ASCII dendrogram:

```sh
python scripts/visualize_dendrogram.py \
  results/dendrograms/average-link/ppop_DATASET_NAME.csv DENDROGRAM_VISUALIZATION_NAME.txt
```

Example:

```sh
python scripts/visualize_dendrogram.py \
  results/dendrograms/average-link/ppop_test_100.csv test_100_vis.txt
```

The output is written to `results/visual/DENDROGRAM_VISUALIZATION_NAME.txt`.

#### Plot cluster colourings by cutting a dendrogram after `n - K_CLUSTERS` merges:

```sh
python scripts/plot_clusters.py --input data/DATASET_NAME.csv --k K_CLUSTERS
python scripts/plot_clusters.py --input data/DATASET_NAME.csv --k K_CLUSTERS --family average --run
```

Example:

```sh
python scripts/plot_clusters.py --input data/test_100.csv --k 4
python scripts/plot_clusters.py --input data/test_100.csv --k 4 --family average --run
```

Useful options:

- `--family single|average|all`
- `--threads N_THREADS`
- `--run`, to generate dendrograms first
- `--out <path>`
- `--dendro name:path`, repeatable

Cluster plots are written to `results/plots/visualizations/clusters_<stem>_k<k>_<family>.png`.

## Benchmarks

#### Run all benchmark families:

```sh
python scripts/benchmark.py
```

#### Run one family:

```sh
python scripts/benchmark.py --family single
python scripts/benchmark.py --family average
```

Benchmark scripts write CSVs to `results/benchmarks/` and figures to `results/plots/analysis/`. They require `pandas` and `matplotlib`.

#### Regenerate benchmark and cluster plots:

```sh
python scripts/benchmark.py
python scripts/plot_clusters.py --input data/DATASET_NAME.csv --k K_CLUSTERS
```

Example:

```sh
python scripts/benchmark.py
python scripts/plot_clusters.py --input data/test_100.csv --k 4
```

Main benchmark outputs:

- `results/benchmarks/single_link_results.csv`
- `results/benchmarks/average_link_results.csv`
- speedup tables and figures under `results/plots/analysis/`