#pragma once
// ─────────────────────────────────────────────────────────────────────────────
// Parameters:
//   data       — N data points (each a vector<double>, any dimension)
//   n_threads  — number of worker threads (includes the calling thread)
//
// Returns the dendrogram as a vector of MergeEvent
// ─────────────────────────────────────────────────────────────────────────────

#include "../sequential/average-link-matrix-seq.hpp"

std::vector<MergeEvent>
hac_average_link_naive(const std::vector<std::vector<double>> &data,
                       int n_threads);
