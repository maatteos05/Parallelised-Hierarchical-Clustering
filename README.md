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

- `data/`: inputs, dendrograms, benchmarks, plots, and visual outputs (see `data/README.md`).
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

All commands below are run from the **repository root**. Inputs are point clouds in `data/inputs/`; dendrogram outputs go in `data/dendrograms/` (see `data/DATA_README.md` for the full layout).

### Generating a dataset

Synthetic clustered 2-D data (`generate_data.py` uses `--n` points, `--k` clusters, optional `--seed`):

```sh
python scripts/generate_data.py --n 100 --k 4 --out data/inputs/test_100.csv
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
  data/inputs/test_100.csv data/dendrograms/single-link/baseline_100.csv

./single-link/hac_single_mst \
  data/inputs/test_100.csv data/dendrograms/single-link/mst_100.csv 4
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
  data/inputs/test_100.csv data/dendrograms/average-link/seq_100.csv

./average-link/hac_naive \
  data/inputs/test_100.csv data/dendrograms/average-link/naive_100.csv 4

./average-link/hac_ppop \
  data/inputs/test_100.csv data/dendrograms/average-link/ppop_100.csv 4
```

Each program prints its wall-clock time to stdout when it finishes.

### Check correctness

**Single-link** -- sequential vs parallel must produce the same merge distances and cluster sizes:

```sh
python scripts/validate_single_link.py \
  --seq data/dendrograms/single-link/baseline_100.csv \
  --par data/dendrograms/single-link/mst_100.csv
```

**Average-link** -- compare any C++ dendrogram against a Python reference implementation on the same input:

```sh
python scripts/validate_average_link.py \
  --data data/inputs/test_100.csv \
  --cpp_out data/dendrograms/average-link/seq_100.csv
```

Use the same command with `naive_100.csv` or `ppop_100.csv` as `--cpp_out` to validate the parallel versions.

### 4. Visualize a dendrogram (optional)

Renders an ASCII tree to `data/visual/`:

```sh
python scripts/visualize_dendrogram.py \
  data/dendrograms/average-link/ppop_100.csv average_ppop_100_vis.txt
```

The output file is written to `data/visual/average_ppop_100_vis.txt`.

### 5. Run performance benchmarks

Runs each implementation over several dataset sizes and thread counts, writes timing CSVs to `data/benchmarks/` and plots to `data/plots/`. Requires `pandas` and `matplotlib`.

```sh
python scripts/single-link-benchmark.py
python scripts/average-link-benchmark.py
```

Results: `data/benchmarks/single_link_results.csv`, `data/benchmarks/average_link_results.csv`, and corresponding speedup tables; figures under `data/plots/`.