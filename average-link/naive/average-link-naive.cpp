#include "average-link-naive.hpp"

#include <atomic>
#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <thread>

static double euclidean(const std::vector<double> &a,
                        const std::vector<double> &b) {
  assert(a.size() == b.size());
  double sum = 0.0;
  for (size_t d = 0; d < a.size(); ++d) {
    double diff = a[d] - b[d];
    sum += diff * diff;
  }
  return std::sqrt(sum);
}

struct LocalMin {
  double dist = std::numeric_limits<double>::infinity();
  int i = -1;
  int j = -1;
};

// ─── Phase 0 — parallel distance-matrix initialisation ──────────────────────
static void init_D_thread(std::atomic<int> &next_row, int N,
                          const std::vector<std::vector<double>> &data,
                          std::vector<std::vector<double>> &D) {
  int i;
  while ((i = next_row.fetch_add(1)) < N) {
    for (int j = i + 1; j < N; ++j) {
      double d = euclidean(data[i], data[j]);
      D[i][j] = D[j][i] = d;
    }
  }
}

// ─── Phase A — parallel find-min ────────────────────────────────────────────
static void find_min_thread(std::atomic<int> &next_row, int n_ids,
                            const std::vector<std::vector<double>> &D,
                            const std::vector<char> &active, LocalMin &out) {
  double best = std::numeric_limits<double>::infinity();
  int best_i = -1;
  int best_j = -1;

  int i;
  while ((i = next_row.fetch_add(1)) < n_ids) {
    if (!active[i])
      continue;
    const std::vector<double> &row = D[i];
    for (int j = i + 1; j < n_ids; ++j) {
      if (!active[j])
        continue;
      double d = row[j];
      if (d < best) {
        best = d;
        best_i = i;
        best_j = j;
      }
    }
  }

  out.dist = best;
  out.i = best_i;
  out.j = best_j;
}

// ─── Phase C — parallel Lance-Williams update ───────────────────────────────
static void update_D_thread(int begin, int end, const std::vector<int> &targets,
                            int k, int i, int j, int sz_i, int sz_j, int sz_k,
                            std::vector<std::vector<double>> &D) {
  for (int idx = begin; idx < end; ++idx) {
    int m = targets[idx];
    double d_km = (sz_i * D[i][m] + sz_j * D[j][m]) / static_cast<double>(sz_k);
    D[k][m] = D[m][k] = d_km;
  }
}

// ─── Main algorithm ─────────────────────────────────────────────────────────
std::vector<MergeEvent>
hac_average_link_naive(const std::vector<std::vector<double>> &data,
                       int n_threads) {

  const int N = static_cast<int>(data.size());
  if (N < 2)
    throw std::invalid_argument("Need at least 2 data points.");
  if (n_threads < 1)
    throw std::invalid_argument("n_threads must be >= 1");

  const int MAX_ID = 2 * N;

  std::vector<std::vector<double>> D(MAX_ID, std::vector<double>(MAX_ID, 0.0));

  // Cluster sizes
  std::vector<int> sz(MAX_ID, 0);
  for (int i = 0; i < N; ++i)
    sz[i] = 1;

  // Active flag per cluster id
  std::vector<char> active(MAX_ID, 0);
  for (int i = 0; i < N; ++i)
    active[i] = 1;

  // ── Phase 0 — parallel D initialisation ──────────────────────────────────
  {
    std::atomic<int> next_row(0);
    std::vector<std::thread> workers;
    workers.reserve(n_threads - 1);
    for (int t = 0; t < n_threads - 1; ++t) {
      workers.emplace_back(init_D_thread, std::ref(next_row), N, std::ref(data),
                           std::ref(D));
    }
    init_D_thread(next_row, N, data, D);
    for (auto &w : workers)
      w.join();
  }

  std::vector<MergeEvent> dendrogram;
  dendrogram.reserve(N - 1);

  int next_id = N;  // id to assign to the next merged cluster
  int n_active = N; // number of currently-active clusters

  // ── Main loop: N-1 merges ───────────────────────────────────────────────
  while (n_active > 1) {

    // ── Phase A — parallel find-min ──────────────────────────────────────
    std::atomic<int> next_row(0);
    std::vector<LocalMin> locals(n_threads);
    std::vector<std::thread> workers;
    workers.reserve(n_threads - 1);

    for (int t = 0; t < n_threads - 1; ++t) {
      workers.emplace_back(find_min_thread, std::ref(next_row), next_id,
                           std::ref(D), std::ref(active), std::ref(locals[t]));
    }
    find_min_thread(next_row, next_id, D, active, locals[n_threads - 1]);
    for (auto &w : workers)
      w.join();

    // Sequential reduction of per-thread mins to the global minimum.
    double best_dist = std::numeric_limits<double>::infinity();
    int best_i = -1;
    int best_j = -1;
    for (const LocalMin &lm : locals) {
      if (lm.dist < best_dist) {
        best_dist = lm.dist;
        best_i = lm.i;
        best_j = lm.j;
      }
    }

    assert(best_i != -1 && best_j != -1);

    // ── Phase B — sequential merge bookkeeping ──────────────────────────
    const int k = next_id++;
    sz[k] = sz[best_i] + sz[best_j];

    std::vector<int> targets;
    targets.reserve(n_active - 2);
    for (int m = 0; m < k; ++m) {
      if (active[m] && m != best_i && m != best_j)
        targets.push_back(m);
    }

    dendrogram.push_back({best_i, best_j, best_dist, sz[k]});

    // ── Phase C — parallel Lance-Williams update ────────────────────────
    const int n_targets = static_cast<int>(targets.size());
    if (n_targets > 0) {
      const int base = n_targets / n_threads;
      const int rem = n_targets % n_threads;

      std::vector<std::thread> upd_workers;
      upd_workers.reserve(n_threads - 1);

      int cursor = 0;
      auto chunk_end = [&](int t) { return cursor + base + (t < rem ? 1 : 0); };

      // Workers handle the first n_threads-1 chunks.
      for (int t = 0; t < n_threads - 1; ++t) {
        int begin = cursor;
        int end = chunk_end(t);
        cursor = end;
        if (begin < end) {
          upd_workers.emplace_back(update_D_thread, begin, end,
                                   std::ref(targets), k, best_i, best_j,
                                   sz[best_i], sz[best_j], sz[k], std::ref(D));
        }
      }
      // Main thread handles the final chunk.
      {
        int begin = cursor;
        int end = n_targets;
        if (begin < end) {
          update_D_thread(begin, end, targets, k, best_i, best_j, sz[best_i],
                          sz[best_j], sz[k], D);
        }
      }
      for (auto &w : upd_workers)
        w.join();
    }

    // Finalise the merge
    active[best_i] = 0;
    active[best_j] = 0;
    active[k] = 1;
    n_active -= 1;
  }

  return dendrogram;
}
