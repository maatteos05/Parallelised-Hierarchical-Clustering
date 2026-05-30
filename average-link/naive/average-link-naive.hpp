#pragma once

#include "../sequential/average-link-matrix-seq.hpp"

std::vector<MergeEvent>
hac_average_link_naive(const std::vector<std::vector<double>> &data,
                       int n_threads);
