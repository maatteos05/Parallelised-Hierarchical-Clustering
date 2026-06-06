#include "average-link-matrix-seq.hpp"

#include <cassert>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>

std::vector<std::vector<double>> load_csv(const std::string &path) {
  std::ifstream file(path);
  if (!file.is_open())
    throw std::runtime_error("Cannot open file: " + path);

  std::vector<std::vector<double>> data;
  std::string line;
  while (std::getline(file, line)) {
    if (line.empty())
      continue;
    std::vector<double> row;
    std::istringstream ss(line);
    std::string token;
    while (std::getline(ss, token, ','))
      row.push_back(std::stod(token));
    data.push_back(std::move(row));
  }
  return data;
}

void save_dendrogram(const std::vector<MergeEvent> &dg, const std::string &path) {
  std::ofstream file(path);
  if (!file.is_open())
    throw std::runtime_error("Cannot open file: " + path);
  file << "cl1,cl2,dist,new_size\n";
  for (const auto &e : dg)
    file << e.cl1 << ',' << e.cl2 << ',' << e.dist << ',' << e.new_size << '\n';
}

static double euclidean(const std::vector<double> &a, const std::vector<double> &b) {
  assert(a.size() == b.size());
  double sum = 0.0;
  for (size_t d = 0; d < a.size(); ++d) {
    double diff = a[d] - b[d];
    sum += diff * diff;
  }
  return std::sqrt(sum);
}

std::vector<MergeEvent>
hac_average_link_matrix(const std::vector<std::vector<double>> &data) {

  const int N = data.size();
  if (N < 2)
    throw std::invalid_argument("Need at least 2 data points.");

  const int MAX_ID = 2 * N;

  std::vector<std::vector<double>> D(MAX_ID, std::vector<double>(MAX_ID, 0.0));

  std::vector<int> sz(MAX_ID, 0);
  for (int i = 0; i < N; ++i)
    sz[i] = 1;

  std::vector<char> active(MAX_ID, 0);
  for (int i = 0; i < N; ++i)
    active[i] = 1;

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

  while (n_active > 1) {

    double best_dist = std::numeric_limits<double>::infinity();
    int best_i = -1;
    int best_j = -1;

    for (int i = 0; i < next_id; ++i) {
      if (!active[i])
        continue;
      for (int j = i + 1; j < next_id; ++j) {
        if (!active[j])
          continue;
        double d = D[i][j];
        if (d < best_dist) {
          best_dist = d;
          best_i = i;
          best_j = j;
        }
      }
    }

    assert(best_i != -1 && best_j != -1);

    const int k = next_id++;
    sz[k] = sz[best_i] + sz[best_j];

    for (int m = 0; m < k; ++m) {
      if (!active[m] || m == best_i || m == best_j)
        continue;
      double d_km =
          (sz[best_i] * D[best_i][m] + sz[best_j] * D[best_j][m]) / (sz[k]);
      D[k][m] = D[m][k] = d_km;
    }

    dendrogram.push_back({best_i, best_j, best_dist, sz[k]});

    active[best_i] = 0;
    active[best_j] = 0;
    active[k] = 1;
    n_active -= 1;
  }

  return dendrogram;
}
