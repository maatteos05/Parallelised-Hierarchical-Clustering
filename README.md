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