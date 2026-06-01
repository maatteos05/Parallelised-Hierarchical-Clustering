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
- `scripts/`: shared Python scripts for data generation (`generate_data.py`) and output validation (`validate_single_link.py`, `validate_average_link.py`).
- `average-link/`: average-link implementations.
- `single-link/`: single-link implementations and notes.

### Build

From the root of the repository, compile every C++ implementation with:

```sh
make
```

This currently builds:

- `average-link/hac_seq`
- `average-link/hac_ppop`
- `single-link/hac_single_baseline`
- `single-link/hac_single_mst`

To build only one part, run:

```sh
make average
make single
```

### Test

Generate a small synthetic dataset:

```sh
python scripts/generate_data.py --n 100 --k 4 --out data/inputs/test_100.csv
```

Run and test single-link:

```sh
./single-link/hac_single_baseline data/inputs/test_100.csv data/dendrograms/single-link/baseline_100.csv
./single-link/hac_single_mst data/inputs/test_100.csv data/dendrograms/single-link/mst_100.csv 4
```

Compare the single-link merge distances and cluster sizes:

```sh
python scripts/validate_single_link.py --seq data/dendrograms/single-link/baseline_100.csv --par data/dendrograms/single-link/mst_100.csv
```

Run and test average-link sequential:

```sh
./average-link/hac_seq data/inputs/test_100.csv data/dendrograms/average-link/seq_100.csv
```

Compare average-link against the Python reference implementation:

```sh
python scripts/validate_average_link.py --data data/inputs/test_100.csv --cpp_out data/dendrograms/average-link/seq_100.csv
```

Run average-link pPOP:

```sh
./average-link/hac_ppop data/inputs/test_100.csv data/dendrograms/average-link/ppop_100.csv [num_threads]
```

Clean compiled binaries with:

```sh
make clean
```