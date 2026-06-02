# Parallelised-Hierarchical-Clustering

## Team: Matteo Sainton, Mateo Fatas, Emeric Payer

### Project scope and decisions

We have decided to implement and parallelize three hierarchical agglomerative clustering (HAC) algorithms:

- Single-link, parallelized via a parallel minimum spanning tree construction, following the PRAM blueprint from Olson (1995) adapted to a CPU multithreaded setting.
- Average-link, parallelized via a thread-parallel distance matrix approach, also following Olson (1995).
- Average-link with pPOP, a second parallelization strategy for average-link using the partially overlapping partitioning scheme from Dash, Petrutiu & Scheuermann (2007). This will allow a direct comparison between a naive parallel approach and a partition-based one on the same linkage criterion.


As a sequential baseline, we will implement the standard O(n²) HAC algorithm for both linkage criteria. All implementations will be in C++ using std::thread primitives consistent with the course material. Benchmarks will be run on datasets of varying sizes (synthetic data generated with controlled cluster structure, as well as datasets from the repository provided in the project description), measuring speedup as a function of thread count.

### Tentative task partition

- Emeric Payer: sequential baseline + parallel single-link (MST-based approach)
- Mateo Fatas: parallel average-link (naive thread-parallel distance matrix)
- Matteo Sainton: pPOP implementation on top of average-link + benchmarking infrastructure

### Repository layout

- `data/`: input point-cloud CSVs (see `data/DATA_README.md`).
- `results/`: dendrograms, benchmark tables, plots, and visual outputs (see `results/README.md`).
- `scripts/`: data generation, validation, dendrogram visualization, and benchmark drivers (see **How to run** below).
- `average-link/`: average-link implementations.
- `single-link/`: single-link implementations and notes.

### Build

From the root of the repository, compile every C++ implementation with:

```sh
make
```

This currently builds:

- `average-link/hac_seq` -- sequential average-link baseline
- `average-link/hac_naive` -- naive parallel average-link
- `average-link/hac_ppop` -- pPOP parallel average-link
- `single-link/hac_single_baseline` -- sequential single-link baseline
- `single-link/hac_single_mst` -- parallel single-link (MST / Borůvka)

To build only one part:

```sh
make average   # average-link only
make single    # single-link only
make clean     # remove all binaries
```

---

## How to run

All commands below are run from the **repository root**. Inputs live in `data/`; outputs go under `results/` (see `data/DATA_README.md` and `results/README.md`).

### Generating a dataset

Synthetic clustered 2-D data (`generate_data.py` uses `--n` points, `--k` clusters, optional `--seed`):

```sh
python scripts/generate_data.py --n 100 --k 4 --out data/test_100.csv
```

Pre-generated inputs already in the repo include `test_100.csv` and `bench_{500,1000,2000,5000}.csv`.

### Run an algorithm

Each binary reads an input CSV and writes a dendrogram CSV (`cl1,cl2,dist,new_size`). Replace paths and thread counts as needed.

**Single-link**

```sh
# Sequential baseline
./single-link/hac_single_baseline <input.csv> <output.csv>

# Parallel (MST) (last argument is thread count)
./single-link/hac_single_mst <input.csv> <output.csv> <threads>
```

Example on the small test set:

```sh
./single-link/hac_single_baseline \
  data/test_100.csv results/dendrograms/single-link/baseline_100.csv

./single-link/hac_single_mst \
  data/test_100.csv results/dendrograms/single-link/mst_100.csv 4
```

**Average-link**

```sh
# Sequential baseline
./average-link/hac_seq <input.csv> <output.csv>

# Naive parallel (optional third argument is thread count (default 4))
./average-link/hac_naive <input.csv> <output.csv> [threads]

# pPOP parallel (third argument is thread count; optional fourth is grid cells per dim)
./average-link/hac_ppop <input.csv> <output.csv> <threads> [n_cells_per_dim]
```

Example:

```sh
./average-link/hac_seq \
  data/test_100.csv results/dendrograms/average-link/seq_100.csv

./average-link/hac_naive \
  data/test_100.csv results/dendrograms/average-link/naive_100.csv 4

./average-link/hac_ppop \
  data/test_100.csv results/dendrograms/average-link/ppop_100.csv 4
```

Each program prints its wall-clock time to stdout when it finishes.

### Check correctness

Runs all binaries (unless `--no-run`), then checks single-link baseline vs MST and each average-link output against a Python reference:

```sh
python scripts/validate.py --input data/test_100.csv
```

Options: `--family single|average|all`, `--threads 4`, `--no-run` (use existing dendrograms).

### Visualize a dendrogram (optional)

Renders an ASCII tree to `results/visual/`:

```sh
python scripts/visualize_dendrogram.py \
  results/dendrograms/average-link/ppop_100.csv average_ppop_100_vis.txt
```

The output file is written to `results/visual/average_ppop_100_vis.txt`.

### Run performance benchmarks

Runs implementations over several dataset sizes and thread counts; writes CSVs to `results/benchmarks/` and plots to `results/plots/`. Requires `pandas` and `matplotlib`.

```sh
python scripts/benchmark.py              # single-link + average-link
python scripts/benchmark.py --family single
python scripts/benchmark.py --family average
```

Results: `results/benchmarks/single_link_results.csv`, `results/benchmarks/average_link_results.csv`, speedup tables, and figures under `results/plots/`.