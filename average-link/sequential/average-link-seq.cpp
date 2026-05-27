#include "average-link-seq.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <limits>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

// ─── I/O ─────────────────────────────────────────────────────────────────────

std::vector<std::vector<double>> load_csv(const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open())
        throw std::runtime_error("Cannot open file: " + path);

    std::vector<std::vector<double>> data;
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::vector<double> row;
        std::istringstream ss(line);
        std::string token;
        while (std::getline(ss, token, ','))
            row.push_back(std::stod(token));
        data.push_back(std::move(row));
    }
    return data;
}

void save_dendrogram(const std::vector<MergeEvent>& dg,
                     const std::string& path) {
    std::ofstream file(path);
    if (!file.is_open())
        throw std::runtime_error("Cannot open file: " + path);
    file << "cl1,cl2,dist,new_size\n";
    for (const auto& e : dg)
        file << e.cl1 << ',' << e.cl2 << ',' << e.dist << ',' << e.new_size << '\n';
}

// ─── Internals ───────────────────────────────────────────────────────────────

// Euclidean distance between two points.
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

// Entry stored inside a per-cluster min-heap.
// Represents: "my distance to cluster `other` is `dist`."
struct PQEntry {
    double dist;
    int    other;
    // min-heap: smallest dist at top
    bool operator>(const PQEntry& o) const { return dist > o.dist; }
};

using MinHeap = std::priority_queue<PQEntry,
                                    std::vector<PQEntry>,
                                    std::greater<PQEntry>>;

// ─── Algorithm ───────────────────────────────────────────────────────────────

std::vector<MergeEvent> hac_average_link(const std::vector<std::vector<double>>& data) {

    const int N = (data.size());
    if (N < 2)
        throw std::invalid_argument("Need at least 2 data points.");

    const int MAX_ID = 2 * N;

    std::vector<std::vector<double>> D(MAX_ID,std::vector<double>(MAX_ID, 0.0));

    std::vector<int> sz(MAX_ID, 0);
    for (int i = 0; i < N; ++i) sz[i] = 1;

    std::vector<MinHeap> pq(MAX_ID);

    // ── Active set ────────────────────────────────────────────────────────
    std::unordered_set<int> active;
    active.reserve(N);
    for (int i = 0; i < N; ++i) active.insert(i);

    // ── Initialise D and heaps ────────────────────────────────────────────
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            double d = euclidean(data[i], data[j]);
            D[i][j] = D[j][i] = d;
            pq[i].push({d, j});
            pq[j].push({d, i});
        }
    }

    // ── Main loop: N-1 merges ─────────────────────────────────────────────
    std::vector<MergeEvent> dendrogram;
    dendrogram.reserve(N - 1);

    int next_id = N;   // id to assign to the next merged cluster

    while (active.size() > 1) {

        double best_dist = std::numeric_limits<double>::infinity();
        int    best_i    = -1;
        int    best_j    = -1;

        for (int i : active) {
            // Discard stale entries (clusters that have been merged away).
            while (!pq[i].empty() && !active.count(pq[i].top().other))
                pq[i].pop();

            if (pq[i].empty()) continue;

            const PQEntry& top = pq[i].top();
            if (top.dist < best_dist) {
                best_dist = top.dist;
                best_i    = i;
                best_j    = top.other;
            }
        }

        assert(best_i != -1 && best_j != -1);

        const int k = next_id++;
        sz[k] = sz[best_i] + sz[best_j];

        for (int m : active) {
            if (m == best_i || m == best_j) continue;

            // Lance-Williams for average-link:
            //   D[k][m] = (|i|·D[i][m] + |j|·D[j][m]) / (|i| + |j|)
            double d_km = ((sz[best_i]) * D[best_i][m] + (sz[best_j]) * D[best_j][m]) / (sz[k]);

            D[k][m] = D[m][k] = d_km;
            pq[k].push({d_km, m});
            pq[m].push({d_km, k});
        }

        dendrogram.push_back({best_i, best_j, best_dist, sz[k]});

        active.erase(best_i);
        active.erase(best_j);
        active.insert(k);
    }

    return dendrogram;
}