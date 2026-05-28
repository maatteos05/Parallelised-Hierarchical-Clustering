// ─────────────────────────────────────────────────────────────────────────────
// Matrix-based sequential baseline for average-link HAC.
//
// See average-link-matrix-seq.hpp for a full rationale of what was changed
// vs. Matteo's heap-based version and why. Short version:
//
//   - Both parallel implementations we will benchmark (naive thread-parallel
//     distance matrix, and pPOP) operate on the distance matrix directly.
//   - To get a clean speedup measurement, the sequential baseline must do
//     the same per-step work on the same data structure — otherwise the
//     reported speedup conflates an algorithmic change (heap -> matrix scan)
//     with the actual parallelization gain.
//   - This file therefore reuses Matteo's I/O helpers, MergeEvent struct,
//     Euclidean distance and Lance-Williams update, but replaces the
//     per-cluster priority queues + unordered_set<active> with a plain
//     matrix scan over a contiguous `active` flag vector.
// ─────────────────────────────────────────────────────────────────────────────

#include "average-link-matrix-seq.hpp"

#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>

// Euclidean distance between two points. Same as in the heap version; kept
// local (static) so the two translation units don't clash at link time.
static double euclidean(const std::vector<double>& a,
                        const std::vector<double>& b) {
    assert(a.size() == b.size());
    double sum = 0.0;
    for (size_t d = 0; d < a.size(); ++d) {
        double diff = a[d] - b[d];
        sum += diff * diff;
    }
    return std::sqrt(sum);
}

std::vector<MergeEvent> hac_average_link_matrix(
    const std::vector<std::vector<double>>& data) {

    const int N = static_cast<int>(data.size());
    if (N < 2)
        throw std::invalid_argument("Need at least 2 data points.");

    // Cluster ids: leaves 0..N-1, merged clusters N..2N-2.
    const int MAX_ID = 2 * N;

    // Distance matrix indexed by cluster id (symmetric). Rows/cols for
    // not-yet-created clusters are written before they are read, so the
    // initial zero-fill is fine.
    std::vector<std::vector<double>> D(MAX_ID, std::vector<double>(MAX_ID, 0.0));

    // Cluster sizes (number of original points in the cluster).
    std::vector<int> sz(MAX_ID, 0);
    for (int i = 0; i < N; ++i) sz[i] = 1;

    // Active flag per cluster id. Using a flat vector<char> (rather than
    // unordered_set<int>) so the find-min scan is a tight, cache-friendly
    // linear pass — same layout the parallel versions will use.
    std::vector<char> active(MAX_ID, 0);
    for (int i = 0; i < N; ++i) active[i] = 1;

    // Initialise upper triangle of D (and mirror to lower triangle).
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            double d = euclidean(data[i], data[j]);
            D[i][j] = D[j][i] = d;
        }
    }

    std::vector<MergeEvent> dendrogram;
    dendrogram.reserve(N - 1);

    int next_id = N;
    int n_active = N;

    // Main loop: N-1 merges.
    while (n_active > 1) {

        // ── Find globally closest pair by scanning the active submatrix ───
        double best_dist = std::numeric_limits<double>::infinity();
        int    best_i    = -1;
        int    best_j    = -1;

        for (int i = 0; i < next_id; ++i) {
            if (!active[i]) continue;
            for (int j = i + 1; j < next_id; ++j) {
                if (!active[j]) continue;
                double d = D[i][j];
                if (d < best_dist) {
                    best_dist = d;
                    best_i    = i;
                    best_j    = j;
                }
            }
        }

        assert(best_i != -1 && best_j != -1);

        // ── Merge best_i and best_j into a new cluster k ──────────────────
        const int k = next_id++;
        sz[k] = sz[best_i] + sz[best_j];

        // Lance-Williams update for average linkage:
        //   D[k][m] = (|i|·D[i][m] + |j|·D[j][m]) / (|i| + |j|)
        for (int m = 0; m < k; ++m) {
            if (!active[m] || m == best_i || m == best_j) continue;
            double d_km = (sz[best_i] * D[best_i][m] + sz[best_j] * D[best_j][m])
                        / static_cast<double>(sz[k]);
            D[k][m] = D[m][k] = d_km;
        }

        dendrogram.push_back({best_i, best_j, best_dist, sz[k]});

        active[best_i] = 0;
        active[best_j] = 0;
        active[k]      = 1;
        n_active -= 1; // removed 2, added 1
    }

    return dendrogram;
}
