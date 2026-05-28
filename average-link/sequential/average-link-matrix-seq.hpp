#pragma once
// ─────────────────────────────────────────────────────────────────────────────
// Matrix-based sequential baseline for average-link HAC.
//
// What changed vs. average-link-seq.{hpp,cpp} (Matteo's heap version) and why:
//
//   1. Find-min strategy: the heap version maintained one std::priority_queue
//      per cluster with lazy deletion of stale entries. This baseline instead
//      finds the closest pair by a plain double loop over the active region
//      of the distance matrix D.
//
//      Why: the two parallel implementations we are benchmarking (naive
//      thread-parallel distance matrix, and pPOP) both operate directly on
//      the distance matrix. Their parallel speedups should be measured
//      against a sequential algorithm with the *same data structure and
//      same per-step work*, otherwise the reported speedup mixes
//      "algorithmic change (heap -> matrix scan)" with the actual
//      parallelization gain. Using a matrix-based baseline gives a clean
//      apples-to-apples comparison.
//
//   2. Removed per-cluster priority queues entirely (no MinHeap, no PQEntry,
//      no lazy-deletion bookkeeping). The matrix D and an `active` flag per
//      cluster are sufficient.
//
//   3. `active` is now a std::vector<char> indexed by cluster id rather than
//      an std::unordered_set<int>. This matches what the parallel versions
//      will use (cheap O(1) check, cache-friendly linear scan) and avoids
//      hash-set overhead inside the inner loops.
//
// Reused unchanged from Matteo's file (via the shared header types):
//   - MergeEvent struct
//   - load_csv / save_dendrogram I/O helpers
//   - Euclidean distance
//   - Lance-Williams update formula for average linkage
// ─────────────────────────────────────────────────────────────────────────────

#include "average-link-seq.hpp"

// Matrix-based sequential average-link HAC.
// Same input/output contract as hac_average_link in average-link-seq.hpp.
std::vector<MergeEvent> hac_average_link_matrix(
    const std::vector<std::vector<double>>& data);
