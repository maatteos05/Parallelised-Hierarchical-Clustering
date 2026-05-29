#include "average-link-pPOP.hpp"

#include <atomic>
#include <barrier>
#include <cmath>
#include <limits>
#include <mutex>
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

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

// True if point (px, py) falls inside the cell boundary extended by delta.
static bool in_cell(double px, double py, const Cell &cell, double delta) {
  return px >= cell.x_min - delta && px <= cell.x_max + delta &&
         py >= cell.y_min - delta && py <= cell.y_max + delta;
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

// Remove value val from vector v in-place.
static void remove_from(std::vector<int> &v, int val) {
  std::vector<int> tmp;
  for (int x : v)
    if (x != val)
      tmp.push_back(x);
  v = tmp;
}

// True if val is already in vector v.
static bool contains(const std::vector<int> &v, int val) {
  for (int x : v)
    if (x == val)
      return true;
  return false;
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

  // Build axis-parallel grid ──────────────────────────────────────────
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
      int cid = row * n_cells_per_dim + col;
      cells[cid].x_min = x_min + col * cell_w;
      cells[cid].x_max = x_min + (col + 1) * cell_w;
      cells[cid].y_min = y_min + row * cell_h;
      cells[cid].y_max = y_min + (row + 1) * cell_h;
    }
  }

  // Assign original points to cells ───────────────────────────────────
  for (int i = 0; i < N; i++) {
    for (int c = 0; c < total_cells; c++) {
      if (in_cell(data[i][0], data[i][1], cells[c], delta)) {
        cells[c].cluster_ids.push_back(i);
        cluster_cells[i].push_back(c);
      }
    }
  }

  // Shared state for the parallel Phase 1 loop ────────────────────────
  double global_min_dist = std::numeric_limits<double>::infinity();
  int global_min_i = -1;
  int global_min_j = -1;
  std::mutex global_min_mutex;

  std::atomic<int> next_cell(0);

  std::atomic<bool> phase1_done(false);

  std::barrier barrier_find_min(n_threads);
  std::barrier barrier_merge(n_threads);

  std::vector<MergeEvent> dendrogram;
  int next_id = N;

  auto worker = [&](int thread_id) {
    while (!phase1_done.load()) {

      double local_min = std::numeric_limits<double>::infinity();
      int local_i = -1;
      int local_j = -1;

      int c;
      while ((c = next_cell.fetch_add(1)) < total_cells) {
        const std::vector<int> &ids = cells[c].cluster_ids;
        for (int a = 0; a < (int)ids.size(); a++) {
          int ci = ids[a];
          if (!active[ci])
            continue;
          for (int b = a + 1; b < (int)ids.size(); b++) {
            int cj = ids[b];
            if (!active[cj])
              continue;
            if (D[ci][cj] < local_min) {
              local_min = D[ci][cj];
              local_i = ci;
              local_j = cj;
            }
          }
        }
      }

      {
        std::lock_guard<std::mutex> lock(global_min_mutex);
        if (local_min < global_min_dist) {
          global_min_dist = local_min;
          global_min_i = local_i;
          global_min_j = local_j;
        }
      }

      barrier_find_min.arrive_and_wait();

      if (thread_id == 0) {

        if (global_min_i == -1 || global_min_dist >= delta) {
          phase1_done.store(true);

        } else {
          int bi = global_min_i;
          int bj = global_min_j;

          // Create new cluster k.
          int k = next_id++;
          sz[k] = sz[bi] + sz[bj];

          std::vector<int> k_cells;
          for (int cell_id : cluster_cells[bi])
            if (!contains(k_cells, cell_id))
              k_cells.push_back(cell_id);
          for (int cell_id : cluster_cells[bj])
            if (!contains(k_cells, cell_id))
              k_cells.push_back(cell_id);
          cluster_cells[k] = k_cells;

          // Lance-Williams update for ALL active clusters.
          for (int m = 0; m < next_id; m++) {
            if (!active[m] || m == bi || m == bj)
              continue;
            double d_km =
                ((double)sz[bi] * D[bi][m] + (double)sz[bj] * D[bj][m]) /
                (double)sz[k];
            D[k][m] = D[m][k] = d_km;
          }

          // Update cell cluster lists.
          for (int cell_id : k_cells) {
            remove_from(cells[cell_id].cluster_ids, bi);
            remove_from(cells[cell_id].cluster_ids, bj);
            if (!contains(cells[cell_id].cluster_ids, k))
              cells[cell_id].cluster_ids.push_back(k);
          }

          // Record merge.
          MergeEvent event;
          event.cl1 = bi;
          event.cl2 = bj;
          event.dist = global_min_dist;
          event.new_size = sz[k];
          dendrogram.push_back(event);

          active[bi] = 0;
          active[bj] = 0;
          active[k] = 1;

          // Reset for next iteration.
          global_min_dist = std::numeric_limits<double>::infinity();
          global_min_i = -1;
          global_min_j = -1;
          next_cell.store(0);
        }
      }

      barrier_merge.arrive_and_wait();
    }
  };

  // Launch threads ────────────────────────────────────────────────────
  std::vector<std::thread> threads;
  for (int t = 0; t < n_threads; t++)
    threads.push_back(std::thread(worker, t));
  for (auto &t : threads)
    t.join();

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