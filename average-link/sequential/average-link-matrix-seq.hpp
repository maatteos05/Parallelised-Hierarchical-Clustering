#pragma once
#include <string>
#include <vector>

struct MergeEvent {
  int cl1, cl2;
  double dist;
  int new_size;
};

std::vector<std::vector<double>> load_csv(const std::string &path);

void save_dendrogram(const std::vector<MergeEvent> &dg,
                     const std::string &path);

std::vector<MergeEvent>
hac_average_link_matrix(const std::vector<std::vector<double>> &data);
