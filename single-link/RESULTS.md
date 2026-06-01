# Single-Link Results

Correctness: all MST outputs matched the sequential baseline.

Speedup formula:

```text
speedup = baseline_time / mst_4_threads_time
```

## Runtime Results


| n    | seed | baseline (ms) | MST 4 threads (ms) | speedup |
| ---- | ---- | ------------- | ------------------ | ------- |
| 100  | 1    | 0.241         | 0.648              | 0.37x   |
| 100  | 2    | 0.369         | 0.972              | 0.38x   |
| 100  | 3    | 0.349         | 0.734              | 0.48x   |
| 100  | 4    | 0.362         | 0.720              | 0.50x   |
| 100  | 5    | 0.236         | 0.805              | 0.29x   |
| 500  | 1    | 7.827         | 5.660              | 1.38x   |
| 500  | 2    | 6.505         | 5.233              | 1.24x   |
| 500  | 3    | 6.927         | 5.448              | 1.27x   |
| 500  | 4    | 6.439         | 4.740              | 1.36x   |
| 500  | 5    | 8.723         | 5.113              | 1.71x   |
| 1000 | 1    | 29.153        | 19.043             | 1.53x   |
| 1000 | 2    | 37.156        | 23.281             | 1.60x   |
| 1000 | 3    | 38.086        | 28.730             | 1.33x   |
| 1000 | 4    | 39.612        | 24.355             | 1.63x   |
| 1000 | 5    | 33.365        | 18.216             | 1.83x   |
| 2000 | 1    | 168.397       | 104.615            | 1.61x   |
| 2000 | 2    | 170.629       | 107.922            | 1.58x   |
| 2000 | 3    | 179.437       | 104.659            | 1.71x   |
| 2000 | 4    | 138.731       | 101.477            | 1.37x   |
| 2000 | 5    | 128.989       | 87.455             | 1.47x   |


## Averages


| n    | avg baseline (ms) | avg MST 4 threads (ms) | avg speedup |
| ---- | ----------------- | ---------------------- | ----------- |
| 100  | 0.311             | 0.776                  | 0.40x       |
| 500  | 7.284             | 5.239                  | 1.39x       |
| 1000 | 35.474            | 22.725                 | 1.56x       |
| 2000 | 157.237           | 101.226                | 1.55x       |


## Summary

For `n=100`, the MST version is slower because thread overhead dominates.

From `n=500` onward, the 4-thread MST version is consistently faster than the sequential baseline.

Observed average speedup for larger datasets: about `1.4x` to `1.6x`.