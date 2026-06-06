#include "average-link-pPOP.hpp"

#include <atomic>
#include <cmath>
#include <condition_variable>
#include <functional>
#include <limits>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <thread>
#include <vector>

struct Cell {
  double x_min, x_max;
  double y_min, y_max;
  std::vector<int> cluster_ids;
};

struct LocalMin {
  double dist = std::numeric_limits<double>::infinity();
  int ci = -1;
  int cj = -1;
};

struct PQEntry {
  double dist;
  int other;
  bool operator>(const PQEntry &o) const { return dist > o.dist; }
};

using MinHeap =
    std::priority_queue<PQEntry, std::vector<PQEntry>, std::greater<PQEntry>>;

class ThreadPool {
public:
  explicit ThreadPool(int n) : n_threads(n) {
    for (int t = 0; t < n - 1; t++)
      workers.emplace_back([this, t] { worker_loop(t); });
  }

  ~ThreadPool() {
    {
      std::lock_guard<std::mutex> lk(mtx);
      stop = true;
      generation++;
    }
    cv_start.notify_all();
    for (auto &w : workers)
      w.join();
  }

  void run(const std::function<void(int)> &fn) {
    job = &fn;
    {
      std::lock_guard<std::mutex> lk(mtx);
      done = 0;
      generation++;
    }
    cv_start.notify_all();

    (*job)(n_threads - 1);

    std::unique_lock<std::mutex> lk(mtx);
    cv_done.wait(lk, [this] { return done == n_threads - 1; });
  }

private:
  void worker_loop(int tid) {
    int my_gen = 0;
    while (true) {
      std::unique_lock<std::mutex> lk(mtx);
      cv_start.wait(lk, [this, &my_gen] { return generation != my_gen || stop; });
      if (stop)
        return;
      my_gen = generation;
      lk.unlock();

      (*job)(tid);

      lk.lock();
      if (++done == n_threads - 1)
        cv_done.notify_one();
    }
  }

  int n_threads;
  std::vector<std::thread> workers;
  std::mutex mtx;
  std::condition_variable cv_start, cv_done;
  const std::function<void(int)> *job = nullptr;
  int generation = 0;
  int done = 0;
  bool stop = false;
};

static bool in_cell(double px, double py, const Cell &cell) {
  return px >= cell.x_min && px <= cell.x_max && py >= cell.y_min &&
         py <= cell.y_max;
}

static double euclidean(const std::vector<double> &a, const std::vector<double> &b) {
  double sum = 0.0;
  for (int d = 0; d < (int)a.size(); d++) {
    double diff = a[d] - b[d];
    sum += diff * diff;
  }
  return std::sqrt(sum);
}

static bool contains(const std::vector<int> &v, int val) {
  for (int x : v)
    if (x == val)
      return true;
  return false;
}

static std::vector<Cell>
build_grid(int n_cells_per_dim, double delta, double x_min, double x_max, double y_min, double y_max, const std::vector<std::vector<double>> &data, int N, int next_id, const std::vector<char> &active, std::vector<std::vector<int>> &cluster_cells, std::vector<double> &rep_x, std::vector<double> &rep_y) {
  int total_cells = n_cells_per_dim * n_cells_per_dim;
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

  for (int i = 0; i < next_id; i++)
    if (active[i])
      cluster_cells[i].clear();

  for (int i = 0; i < next_id; i++) {
    if (!active[i])
      continue;
    double px = (i < N) ? data[i][0] : rep_x[i];
    double py = (i < N) ? data[i][1] : rep_y[i];
    for (int c = 0; c < total_cells; c++) {
      if (in_cell(px, py, cells[c])) {
        cells[c].cluster_ids.push_back(i);
        cluster_cells[i].push_back(c);
      }
    }
  }

  return cells;
}

static void
build_heaps_thread(std::atomic<int> &next_cell, int total_cells, const std::vector<Cell> &cells, const std::vector<std::vector<int>> &cluster_cells, const std::vector<std::vector<double>> &D, std::vector<MinHeap> &pq) {
  int c;
  while ((c = next_cell.fetch_add(1)) < total_cells) {
    const std::vector<int> &ids = cells[c].cluster_ids;
    for (int i : ids) {
      if (cluster_cells[i][0] != c)
        continue;
      for (int cell_id : cluster_cells[i]) {
        for (int j : cells[cell_id].cluster_ids) {
          if (j != i)
            pq[i].push({D[i][j], j});
        }
      }
    }
  }
}

static void find_min_thread(std::atomic<int> &next_cell, int total_cells, const std::vector<Cell> &cells, const std::vector<char> &active, const std::vector<std::vector<int>> &cluster_cells, std::vector<MinHeap> &pq, LocalMin &result) {
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

static void update_heaps_thread(int begin, int end, const std::vector<int> &affected_clusters, int k, const std::vector<std::vector<double>> &D, std::vector<MinHeap> &pq) {
  for (int i = begin; i < end; i++) {
    int m = affected_clusters[i];
    pq[m].push({D[k][m], k});
  }
}

static void update_D_thread(int begin, int end, const std::vector<int> &targets, int k, int best_i, int best_j, int sz_i, int sz_j, int sz_k, std::vector<std::vector<double>> &D) {
  for (int idx = begin; idx < end; idx++) {
    int m = targets[idx];
    double d_km = ((double)sz_i * D[best_i][m] + (double)sz_j * D[best_j][m]) / (double)sz_k;
    D[k][m] = D[m][k] = d_km;
  }
}

static inline int chunk_begin(int tid, int n, int n_threads) {
  int block = n / n_threads;
  return tid * block;
}
static inline int chunk_end(int tid, int n, int n_threads) {
  int block = n / n_threads;
  return (tid == n_threads - 1) ? n : (tid + 1) * block;
}

std::vector<MergeEvent>
hac_average_link_ppop(const std::vector<std::vector<double>> &data, int n_threads, int initial_n_cells_per_dim) {
  int N = (int)data.size();
  int MAX_ID = 2 * N;

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
  int n_active = N;

  std::vector<std::vector<int>> cluster_cells(MAX_ID);

  for (int i = 0; i < N; i++)
    for (int j = i + 1; j < N; j++)
      D[i][j] = D[j][i] = euclidean(data[i], data[j]);

  std::vector<MinHeap> pq(MAX_ID);

  std::vector<double> rep_x(MAX_ID, 0.0);
  std::vector<double> rep_y(MAX_ID, 0.0);
  for (int i = 0; i < N; i++) {
    rep_x[i] = data[i][0];
    rep_y[i] = data[i][1];
  }

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

  const double D_INCR = 0.1;
  const double C_DECR = 0.1;

  double delta = 0.0;
  int n_cells_per_dim = initial_n_cells_per_dim;
  double last_exceeding_dist = 0.0;

  std::vector<MergeEvent> dendrogram;
  int next_id = N;

  std::atomic<int> next_cell(0);

  ThreadPool pool(n_threads);

  while (n_active > 1) {

    int total_cells = n_cells_per_dim * n_cells_per_dim;

    std::vector<Cell> cells =
        build_grid(n_cells_per_dim, delta, x_min, x_max, y_min, y_max, data, N,
                   next_id, active, cluster_cells, rep_x, rep_y);

    for (int i = 0; i < next_id; i++) {
      if (active[i]) {
        MinHeap empty;
        std::swap(pq[i], empty);
      }
    }

    next_cell.store(0);
    pool.run([&](int) {
      build_heaps_thread(next_cell, total_cells, cells, cluster_cells, D, pq);
    });

    while (n_active > 1) {

      next_cell.store(0);

      std::vector<LocalMin> results(n_threads);
      pool.run([&](int tid) {
        find_min_thread(next_cell, total_cells, cells, active, cluster_cells,
                        pq, results[tid]);
      });

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

      if (best_i == -1 || best_dist >= delta) {
        last_exceeding_dist = best_dist;
        break;
      }
      int k = next_id++;
      sz[k] = sz[best_i] + sz[best_j];

      rep_x[k] = ((double)sz[best_i] * rep_x[best_i] +
                  (double)sz[best_j] * rep_x[best_j]) /
                 sz[k];
      rep_y[k] = ((double)sz[best_i] * rep_y[best_i] +
                  (double)sz[best_j] * rep_y[best_j]) /
                 sz[k];

      std::vector<int> k_cells;
      for (int cell_id : cluster_cells[best_i])
        if (!contains(k_cells, cell_id))
          k_cells.push_back(cell_id);
      for (int cell_id : cluster_cells[best_j])
        if (!contains(k_cells, cell_id))
          k_cells.push_back(cell_id);
      cluster_cells[k] = k_cells;

      std::vector<int> targets;
      for (int m = 0; m < next_id; m++)
        if (active[m] && m != best_i && m != best_j)
          targets.push_back(m);

      int n_targets = (int)targets.size();
      pool.run([&](int tid) {
        update_D_thread(chunk_begin(tid, n_targets, n_threads), chunk_end(tid, n_targets, n_threads), targets, k, best_i, best_j, sz[best_i], sz[best_j], sz[k], D);
      });

      for (int cell_id : k_cells) {
        cells[cell_id].cluster_ids.push_back(k);
      }

      active[best_i] = 0;
      active[best_j] = 0;
      active[k] = 1;
      n_active--;

      MergeEvent event;
      event.cl1 = best_i;
      event.cl2 = best_j;
      event.dist = best_dist;
      event.new_size = sz[k];
      dendrogram.push_back(event);

      std::vector<int> affected_clusters;
      for (int cell_id : k_cells) {
        for (int m : cells[cell_id].cluster_ids) {
          if (active[m] && m != k && !contains(affected_clusters, m))
            affected_clusters.push_back(m);
        }
      }

      int n_affected = (int)affected_clusters.size();
      pool.run([&](int tid) {
        update_heaps_thread(chunk_begin(tid, n_affected, n_threads), chunk_end(tid, n_affected, n_threads), affected_clusters, k, D, pq);
      });

      for (int m : affected_clusters) {
        pq[k].push({D[k][m], m});
      }
    }

    if (last_exceeding_dist == std::numeric_limits<double>::infinity())
      break;

    delta = (1.0 + D_INCR) * last_exceeding_dist;

    int decayed = (int)(n_cells_per_dim * (1.0 - C_DECR));
    if (decayed >= n_cells_per_dim)
      decayed = n_cells_per_dim - 1;
    int floor_for_threads = std::max(2, (int)std::sqrt((double)(2 * n_threads)));
    int cap_by_clusters = std::max(1, (int)std::sqrt((double)n_active / 4.0));
    int floor_cells = std::min(floor_for_threads, cap_by_clusters);
    n_cells_per_dim = std::max(decayed, floor_cells);
  }

  return dendrogram;
}