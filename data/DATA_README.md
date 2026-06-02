# Data — input datasets

Point-cloud CSVs fed to HAC binaries: one `x,y` pair per line, no header.

| Pattern | Description |
|---------|-------------|
| `test*.csv` | Small datasets for correctness runs (e.g. `test_100.csv`). |
| `bench_{N}.csv` | Benchmark inputs (N = 500, 1000, 2000, 5000). |
| `bench_n*_k*_s*.csv` | Extra benchmark inputs from earlier local runs. |
| `sl_{N}_{run}.csv` | Single-link multi-seed inputs. |

Regenerate with:

```sh
python scripts/generate_data.py --n 100 --k 4 --out data/test_100.csv
python scripts/generate_data.py --n 500 --k 10 --seed 42 --out data/bench_500.csv
```

Outputs (dendrograms, benchmarks, plots) live under `results/` — see `results/README.md`.
