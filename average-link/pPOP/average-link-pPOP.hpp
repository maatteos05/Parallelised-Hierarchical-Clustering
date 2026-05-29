#pragma once
// ─────────────────────────────────────────────────────────────────────────────
// Parameters:
//   data            — N data points, each a 2-D coordinate vector
//   n_threads       — number of parallel threads to use
//   n_cells_per_dim — grid side length; total cells = n_cells_per_dim^2
//   delta           — overlap margin (same unit as the data coordinates)
//                     Controls the Phase 1 / Phase 2 cutoff:
//                       delta = 0   → Phase 1 does nothing, all work in Phase 2
//                       delta = inf → Phase 1 does everything sequentially
//                     A good starting value: ~10th percentile of pairwise
//                     distances.
// ─────────────────────────────────────────────────────────────────────────────

#include "../sequential/average-link-matrix-seq.hpp"

std::vector<MergeEvent>
hac_average_link_ppop(const std::vector<std::vector<double>> &data,
                      int n_threads, int n_cells_per_dim, double delta);