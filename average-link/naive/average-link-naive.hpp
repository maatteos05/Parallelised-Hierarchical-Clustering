#pragma once
// ─────────────────────────────────────────────────────────────────────────────
// Naive thread-parallel average-link HAC.
//
// Direct parallelization of the matrix-based sequential algorithm (Olson 1995,
// §4.4): same data structures (N×N distance matrix, active flags, sizes),
// same per-step work, just split across threads phase by phase. No heaps,
// no spatial partitioning — this is the baseline-parallel reference against
// which the smarter pPOP implementation is compared.
//
// Per merge step:
//   Phase A — parallel find-min over active pairs
//   Phase B — sequential merge bookkeeping (cheap O(n))
//   Phase C — parallel Lance-Williams update for the new cluster's row/column
//
// Parameters:
//   data       — N data points (each a vector<double>, any dimension)
//   n_threads  — number of worker threads (includes the calling thread)
//
// Returns the dendrogram as a vector of MergeEvent (defined in the shared
// sequential header).
// ─────────────────────────────────────────────────────────────────────────────

#include "../sequential/average-link-matrix-seq.hpp"

std::vector<MergeEvent>
hac_average_link_naive(const std::vector<std::vector<double>> &data,
                       int n_threads);
