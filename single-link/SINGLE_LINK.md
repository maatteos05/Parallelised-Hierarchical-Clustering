# Single-Link

- `baseline_olson_matrix.cpp`: sequential baseline (follows matrix-based idea from Olson rather than recomputing every point distance at every merge)
- `parallel_mst.cpp`: parallel version (uses the fact that single-link HAC can be recovered from a minimum spanning tree, then builds the MST with a threaded Boruvka-style approach)

## Compile

From this folder:

```sh
make
```

Or from the repository root:

```sh
make single
```

This builds:

- `hac_single_baseline`
- `hac_single_mst`

## Run

Baseline:

```sh
./hac_single_baseline ../data/inputs/input.csv ../data/dendrograms/single-link/baseline_output.csv
```

MST version:

```sh
./hac_single_mst ../data/inputs/input.csv ../data/dendrograms/single-link/mst_output.csv 4
```

The last argument is the number of threads.

## Test

Generate data:

```sh
python ../scripts/generate_data.py --n 100 --k 4 --out ../data/inputs/test_100.csv
```

Run both versions:

```sh
./hac_single_baseline ../data/inputs/test_100.csv ../data/dendrograms/single-link/baseline_100.csv
./hac_single_mst ../data/inputs/test_100.csv ../data/dendrograms/single-link/mst_100.csv 4
```

Compare the merge distances and cluster sizes:

```sh
python ../scripts/validate_single_link.py --seq ../data/dendrograms/single-link/baseline_100.csv --par ../data/dendrograms/single-link/mst_100.csv
```
