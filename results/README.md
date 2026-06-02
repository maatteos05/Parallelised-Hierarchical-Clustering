# Results — algorithm outputs

Everything produced by running HAC binaries or analysis scripts.

| Folder | Contents |
|--------|----------|
| `dendrograms/single-link/` | Single-link merge CSVs (`cl1,cl2,dist,new_size`). |
| `dendrograms/average-link/` | Average-link merge CSVs. |
| `benchmarks/` | Timing tables from `scripts/benchmark.py`. |
| `plots/analysis/` | Speedup and runtime PNGs from `scripts/benchmark.py`. |
| `plots/visualizations/` | Cluster colourings from `scripts/plot_clusters.py`. |
| `visual/` | ASCII dendrogram text from `scripts/visualize_dendrogram.py`. |

Regenerate benchmarks and plots:

```sh
python scripts/benchmark.py
python scripts/plot_clusters.py --input data/test_100.csv --k 4
```

## Legacy results

Early naive-only average-link benchmark (different setup, written analysis):

`average-link/naive/results/` — see `BENCHMARK_NOTES.md` and `naive_benchmark.csv`.
