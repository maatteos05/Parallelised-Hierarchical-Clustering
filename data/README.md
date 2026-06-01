# Data layout

| Folder | Contents |
|--------|----------|
| `inputs/` | Point-cloud CSVs (`x,y` per line, no header) fed to HAC binaries. |
| `dendrograms/single-link/` | Single-link merge outputs (`cl1,cl2,dist,new_size`). |
| `dendrograms/average-link/` | Average-link merge outputs. |
| `benchmarks/` | Timing tables from `scripts/*-benchmark.py`. |
| `plots/` | Speedup and runtime PNG figures. |
| `visual/` | ASCII dendrogram text from `scripts/visualize_dendrogram.py`. |

## inputs/

| Pattern | Description |
|---------|-------------|
| `test*.csv` | Small datasets for manual / correctness runs. |
| `bench_{N}.csv` | Benchmark inputs (N = 500, 1000, 2000, 5000). |
| `bench_n*_k*_s*.csv` | Extra benchmark inputs from naive-local runs. |
| `sl_{N}_{run}.csv` | Single-link multi-seed inputs (N = 100, 500, 1000, 2000; run 1–5). |

## dendrograms/single-link/

| Pattern | Description |
|---------|-------------|
| `baseline_*.csv`, `mst_*.csv` | Main single-link outputs. |
| `sl_base_{N}_{run}.csv` | Baseline outputs for multi-seed runs. |
| `sl_mst_{N}_{run}.csv` | MST parallel outputs for multi-seed runs. |

## dendrograms/average-link/

| File | Description |
|------|-------------|
| `seq_100.csv` | Sequential average-link on `test_100`. |
| `ppop_100.csv` | pPOP parallel on `test_100`. |
| `dendrogram_ppop_100.csv` | Alternate pPOP output (distinct from `ppop_100.csv`). |

## benchmarks/

| File | Former name at `data/` root |
|------|----------------------------|
| `single_link_results.csv` | `single_link_benchmark_results.csv` |
| `single_link_speedup.csv` | `single_link_benchmark_speedup.csv` |
| `average_link_results.csv` | `average_link_benchmark_results.csv` |
| `average_link_speedup.csv` | `average_link_benchmark_speedup.csv` |

Benchmark scripts write to these paths. Regenerate with:

```sh
python scripts/single-link-benchmark.py
python scripts/average-link-benchmark.py
```
