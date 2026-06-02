# Results — algorithm outputs

Everything produced by running HAC binaries or analysis scripts.

| Folder | Contents |
|--------|----------|
| `dendrograms/single-link/` | Single-link merge CSVs (`cl1,cl2,dist,new_size`). |
| `dendrograms/average-link/` | Average-link merge CSVs. |
| `benchmarks/` | Timing tables from `scripts/benchmark.py`. |
| `plots/` | Speedup and runtime PNG figures. |
| `visual/` | ASCII dendrogram text from `scripts/visualize_dendrogram.py`. |

Regenerate benchmarks and plots:

```sh
python scripts/benchmark.py
```

## Legacy results

Early naive-only average-link benchmark (different setup, written analysis):

`average-link/naive/results/` — see `BENCHMARK_NOTES.md` and `naive_benchmark.csv`.
