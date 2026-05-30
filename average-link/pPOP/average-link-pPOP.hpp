#pragma once
// ─────────────────────────────────────────────────────────────────────────────
// Parameters:
//   data            — N data points, each a 2-D coordinate vector
//   n_threads       — number of parallel threads to use
//   initial_n_cells_per_dim — grid side length for first iteration;
//                             total initial cells = initial_n_cells_per_dim^2.
//                             Recommended: sqrt(N/10) for ~10 clusters/cell.
// ─────────────────────────────────────────────────────────────────────────────

#include "../sequential/average-link-matrix-seq.hpp"

std::vector<MergeEvent>
hac_average_link_ppop(const std::vector<std::vector<double>> &data,
                      int n_threads, int initial_n_cells_per_dim);