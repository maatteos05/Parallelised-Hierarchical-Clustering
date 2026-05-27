#pragma once
#include <vector>
#include <string>

// ─── Data structures ─────────────────────────────────────────────────────────

// One row of the dendrogram: records a single merge event.
// cl1, cl2  : ids of the two clusters that were merged
// dist      : average-link distance at which they were merged
// new_size  : number of points in the resulting cluster
struct MergeEvent {
    int    cl1, cl2;
    double dist;
    int    new_size;
};

// ─── I/O helpers ─────────────────────────────────────────────────────────────

// Load a CSV file (no header) where each row is a data point.
// Returns a 2-D vector: data[point_index][dimension].
std::vector<std::vector<double>> load_csv(const std::string& path);

// Write the dendrogram to a CSV file (cl1,cl2,dist,new_size per row).
void save_dendrogram(const std::vector<MergeEvent>& dg,
                     const std::string& path);

// ─── Algorithm ───────────────────────────────────────────────────────────────

std::vector<MergeEvent> hac_average_link(const std::vector<std::vector<double>>& data);