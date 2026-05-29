#include "average-link-pPOP.hpp"

#include <atomic>
#include <cmath>
#include <limits>
#include <queue>
#include <stdexcept>
#include <thread>
#include <vector>

// ─────────────────────────────────────────────────────────────────────────────
// Internal data structures
// ─────────────────────────────────────────────────────────────────────────────

struct Cell {
  double x_min, x_max;
  double y_min, y_max;
  std::vector<int> cluster_ids; // ids of clusters currently in this cell
};

// Result of one thread's find-min scan over its assigned cells.
struct LocalMin {
  double dist = std::numeric_limits<double>::infinity();
  int ci = -1;
  int cj = -1;
};

// One entry in a per-cluster min-heap.
// Means: "my distance to cluster `other` is `dist`."
// Entries can become stale (other was merged away) — discarded lazily.
struct PQEntry {
  double dist;
  int other;
  bool operator>(const PQEntry &o) const { return dist > o.dist; }
};

using MinHeap =
    std::priority_queue<PQEntry, std::vector<PQEntry>, std::greater<PQEntry>>;
// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

// True if point (px, py) falls inside the cell boundary extended by delta.
static bool in_cell(double px, double py, const Cell &cell) {
  return px >= cell.x_min && px <= cell.x_max && py >= cell.y_min &&
         py <= cell.y_max;
}

// Euclidean distance between two points (any dimension).
static double euclidean(const std::vector<double> &a,
                        const std::vector<double> &b) {
  double sum = 0.0;
  for (int d = 0; d < (int)a.size(); d++) {
    double diff = a[d] - b[d];
    sum += diff * diff;
  }
  return std::sqrt(sum);
}

// True if val is already in vector v.
static bool contains(const std::vector<int> &v, int val) {
  for (int x : v)
    if (x == val)
      return true;
  return false;
}

// ─────────────────────────────────────────────────────────────────────────────
// Thread worker — scans a fixed chunk of cells [cell_begin, cell_end)
// and finds the closest active pair within those cells.
// ─────────────────────────────────────────────────────────────────────────────
static void find_min_thread(std::atomic<int> &next_cell, int total_cells,
                            const std::vector<Cell> &cells,
                            const std::vector<char> &active,
                            const std::vector<std::vector<int>> &cluster_cells,
                            std::vector<MinHeap> &pq, LocalMin &result) {
  double local_dist = std::numeric_limits<double>::infinity();
  int local_i = -1;
  int local_j = -1;

  int c;
  while ((c = next_cell.fetch_add(1)) < total_cells) {
    const std::vector<int> &ids = cells[c].cluster_ids;
    for (int ci : ids) {
      if (!active[ci])
        continue;
      if (cluster_cells[ci][0] != c)
        continue;
      // Pop stale entries: discard tops that point to inactive clusters.
      // We only check active[] — not cell membership — because a cluster
      // can legitimately belong to multiple cells and its heap entry is
      // valid as long as the other cluster is still active.
      while (!pq[ci].empty() && !active[pq[ci].top().other])
        pq[ci].pop();

      if (pq[ci].empty())
        continue;

      double d = pq[ci].top().dist;
      int cj = pq[ci].top().other;

      if (d < local_dist) {
        local_dist = d;
        local_i = ci;
        local_j = cj;
      }
    }
  }

  result.dist = local_dist;
  result.ci = local_i;
  result.cj = local_j;
}

static void update_heaps_thread(int begin, int end,
                                const std::vector<int> &affected_clusters,
                                int k,
                                const std::vector<std::vector<double>> &D,
                                std::vector<MinHeap> &pq) {
  for (int i = begin; i < end; i++) {
    int m = affected_clusters[i];
    pq[m].push({D[k][m], k});
  }
}
// ─────────────────────────────────────────────────────────────────────────────
// Main algorithm
// ─────────────────────────────────────────────────────────────────────────────

std::vector<MergeEvent>
hac_average_link_ppop(const std::vector<std::vector<double>> &data,
                      int n_threads, int n_cells_per_dim, double delta) {
  int N = (int)data.size();
  int MAX_ID = 2 * N;
  int total_cells = n_cells_per_dim * n_cells_per_dim;

  if (N < 2)
    throw std::invalid_argument("Need at least 2 data points.");
  if ((int)data[0].size() < 2)
    throw std::invalid_argument("pPOP requires at least 2-D data.");

  std::vector<std::vector<double>> D(MAX_ID, std::vector<double>(MAX_ID, 0.0));

  std::vector<int> sz(MAX_ID, 0);
  for (int i = 0; i < N; i++)
    sz[i] = 1;

  std::vector<char> active(MAX_ID, 0);
  for (int i = 0; i < N; i++)
    active[i] = 1;

  std::vector<std::vector<int>> cluster_cells(MAX_ID);

  // Initialise D ──────────────────────────────────────────────────────
  for (int i = 0; i < N; i++)
    for (int j = i + 1; j < N; j++)
      D[i][j] = D[j][i] = euclidean(data[i], data[j]);

  // pq[i] contains (dist, j) entries for every cluster j in the same cell as i.
  std::vector<MinHeap> pq(MAX_ID);

  // Build cells grid ──────────────────────────────────────────
  double x_min = data[0][0], x_max = data[0][0];
  double y_min = data[0][1], y_max = data[0][1];
  for (int i = 0; i < N; i++) {
    if (data[i][0] < x_min)
      x_min = data[i][0];
    if (data[i][0] > x_max)
      x_max = data[i][0];
    if (data[i][1] < y_min)
      y_min = data[i][1];
    if (data[i][1] > y_max)
      y_max = data[i][1];
  }
  x_min -= 1e-9;
  x_max += 1e-9;
  y_min -= 1e-9;
  y_max += 1e-9;

  double cell_w = (x_max - x_min) / n_cells_per_dim;
  double cell_h = (y_max - y_min) / n_cells_per_dim;

  std::vector<Cell> cells(total_cells);
  for (int row = 0; row < n_cells_per_dim; row++) {
    for (int col = 0; col < n_cells_per_dim; col++) {
      int i = row * n_cells_per_dim + col;
      cells[i].x_min = x_min + col * cell_w - delta / 2.0;
      cells[i].x_max = x_min + (col + 1) * cell_w + delta / 2.0;
      cells[i].y_min = y_min + row * cell_h - delta / 2.0;
      cells[i].y_max = y_min + (row + 1) * cell_h + delta / 2.0;
    }
  }

  // Assign original points to cells ─────────────────────────────────────────
  for (int i = 0; i < N; i++) {
    for (int c = 0; c < total_cells; c++) {
      if (in_cell(data[i][0], data[i][1], cells[c])) {
        cells[c].cluster_ids.push_back(i);
        cluster_cells[i].push_back(c);
      }
    }
  }

  // Build heaps: for each cell, for each pair cluster i,j in that cell
  for (int c = 0; c < total_cells; c++) {
    const std::vector<int> &ids = cells[c].cluster_ids;
    for (int a = 0; a < (int)ids.size(); a++) {
      for (int b = a + 1; b < (int)ids.size(); b++) {
        int i = ids[a];
        int j = ids[b];
        pq[i].push({D[i][j], j});
        pq[j].push({D[j][i], i});
      }
    }
  }

  // ─────────────────────────────────────────────────────────────────────
  // Phase 1 — parallel find-min (heap tops) + parallel heap update
  // ─────────────────────────────────────────────────────────────────────

  std::vector<MergeEvent> dendrogram;
  int next_id = N;

  std::atomic<int> next_cell(0);

  while (true) {

    // ── Step 1: parallel find-min via heap tops ───────────────────────
    next_cell.store(0);

    std::vector<LocalMin> results(n_threads);
    std::vector<std::thread> workers(n_threads - 1);

    for (int t = 0; t < n_threads - 1; t++) {
      workers[t] = std::thread(find_min_thread, std::ref(next_cell),
                               total_cells, std::ref(cells), std::ref(active),
                               std::ref(cluster_cells), std::ref(pq),
                               std::ref(results[t]));
    }
    find_min_thread(next_cell, total_cells, cells, active, cluster_cells, pq,
                    results[n_threads - 1]);

    for (int t = 0; t < n_threads - 1; t++)
      workers[t].join();

    // ── Step 2: pick global minimum ───────────────────────────────────
    double best_dist = std::numeric_limits<double>::infinity();
    int best_i = -1;
    int best_j = -1;

    for (int t = 0; t < n_threads; t++) {
      if (results[t].dist < best_dist) {
        best_dist = results[t].dist;
        best_i = results[t].ci;
        best_j = results[t].cj;
      }
    }

    // ── Check termination ─────────────────────────────────────────────
    if (best_i == -1 || best_dist >= delta)
      break;

    // ── Step 3: merge and update D (sequential) ───────────────────────
    int k = next_id++;
    sz[k] = sz[best_i] + sz[best_j];

    // Cell membership of k = union of best_i's and best_j's cells.
    std::vector<int> k_cells;
    for (int cell_id : cluster_cells[best_i])
      if (!contains(k_cells, cell_id))
        k_cells.push_back(cell_id);
    for (int cell_id : cluster_cells[best_j])
      if (!contains(k_cells, cell_id))
        k_cells.push_back(cell_id);
    cluster_cells[k] = k_cells;

    // Lance-Williams update for ALL active clusters.
    for (int m = 0; m < next_id; m++) {
      if (!active[m] || m == best_i || m == best_j)
        continue;
      double d_km = ((double)sz[best_i] * D[best_i][m] +
                     (double)sz[best_j] * D[best_j][m]) /
                    (double)sz[k];
      D[k][m] = D[m][k] = d_km;
    }

    // Update cell cluster lists.
    for (int cell_id : k_cells) {
      cells[cell_id].cluster_ids.push_back(k);
    }

    active[best_i] = 0;
    active[best_j] = 0;
    active[k] = 1;

    // Record merge and update active flags.
    MergeEvent event;
    event.cl1 = best_i;
    event.cl2 = best_j;
    event.dist = best_dist;
    event.new_size = sz[k];
    dendrogram.push_back(event);

    // ── Step 4: parallel heap update ──────────────────────────────────
    // Collect all clusters in the affected cells (deduplication needed
    // since a cluster can appear in multiple cells).
    std::vector<int> affected_clusters;
    for (int cell_id : k_cells) {
      for (int m : cells[cell_id].cluster_ids) {
        if (active[m] && m != k && !contains(affected_clusters, m))
          affected_clusters.push_back(m);
      }
    }

    // Split affected clusters across threads.
    // Each thread pushes (D[k][m], k) into pq[m] for its chunk.
    // No race condition: each thread writes to different pq[m].
    int n_affected = (int)affected_clusters.size();
    int heap_block = (n_affected > 0) ? (n_affected / n_threads) : 0;

    std::vector<std::thread> heap_workers(n_threads - 1);
    int hstart = 0;
    for (int t = 0; t < n_threads - 1; t++) {
      int hend = hstart + heap_block;
      heap_workers[t] = std::thread(update_heaps_thread, hstart, hend,
                                    std::ref(affected_clusters), k, std::ref(D),
                                    std::ref(pq));
      hstart = hend;
    }
    update_heaps_thread(hstart, n_affected, affected_clusters, k, D, pq);

    for (int t = 0; t < n_threads - 1; t++)
      heap_workers[t].join();

    // ── Step 5: build pq[k] from scratch ──────────────────────────────
    // k is brand new — its heap must be populated with distances to
    // all active clusters in its cells.
    for (int m : affected_clusters) {
      pq[k].push({D[k][m], m});
    }
  }

  // Phase 2: sequential HAC on remaining clusters ─────────────────────
  // We finish the dendrogram with a plain matrix-scan over active clusters.
  std::vector<int> remaining;
  for (int i = 0; i < next_id; i++)
    if (active[i])
      remaining.push_back(i);

  int n_remaining = (int)remaining.size();

  while (n_remaining > 1) {

    double best_dist = std::numeric_limits<double>::infinity();
    int best_i = -1;
    int best_j = -1;

    for (int a = 0; a < n_remaining; a++) {
      int i = remaining[a];
      for (int b = a + 1; b < n_remaining; b++) {
        int j = remaining[b];
        if (D[i][j] < best_dist) {
          best_dist = D[i][j];
          best_i = i;
          best_j = j;
        }
      }
    }

    int k = next_id++;
    sz[k] = sz[best_i] + sz[best_j];

    for (int m : remaining) {
      if (m == best_i || m == best_j)
        continue;
      double d_km = ((double)sz[best_i] * D[best_i][m] +
                     (double)sz[best_j] * D[best_j][m]) /
                    (double)sz[k];
      D[k][m] = D[m][k] = d_km;
    }

    MergeEvent event;
    event.cl1 = best_i;
    event.cl2 = best_j;
    event.dist = best_dist;
    event.new_size = sz[k];
    dendrogram.push_back(event);

    active[best_i] = 0;
    active[best_j] = 0;
    active[k] = 1;

    std::vector<int> new_remaining;
    for (int id : remaining)
      if (id != best_i && id != best_j)
        new_remaining.push_back(id);
    new_remaining.push_back(k);
    remaining = new_remaining;
    n_remaining = (int)remaining.size();
  }

  return dendrogram;
}