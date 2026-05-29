#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

struct MergeLine {
    int first_cluster;
    int second_cluster;
    double distance;
    int merged_size;
};

struct WeightedEdge {
    int first_point;
    int second_point;
    double distance;
};

struct Components {
    std::vector<int> parent;
    std::vector<int> tree_height;
    std::vector<int> component_size;

    explicit Components(int number_of_points)
        : parent(number_of_points),
          tree_height(number_of_points, 0),
          component_size(number_of_points, 1) {
        std::iota(parent.begin(), parent.end(), 0);
    }

    int find(int point_id) {
        if (parent[point_id] != point_id) {
            parent[point_id] = find(parent[point_id]);
        }
        return parent[point_id];
    }

    int find_without_compression(int point_id) const {
        while (parent[point_id] != point_id) {
            point_id = parent[point_id];
        }
        return point_id;
    }

    int merge_roots(int first_root, int second_root) {
        if (tree_height[first_root] < tree_height[second_root]) {
            std::swap(first_root, second_root);
        }

        parent[second_root] = first_root;
        component_size[first_root] += component_size[second_root];
        if (tree_height[first_root] == tree_height[second_root]) {
            ++tree_height[first_root];
        }
        return first_root;
    }
};

std::vector<std::vector<double>> read_points_csv(const std::string& file_name) {
    std::ifstream input_file(file_name);
    if (!input_file.is_open()) {
        throw std::runtime_error("Could not open " + file_name);
    }

    std::vector<std::vector<double>> points;
    std::string line;
    while (std::getline(input_file, line)) {
        if (line.empty()) {
            continue;
        }

        std::vector<double> point;
        std::stringstream line_stream(line);
        std::string value;
        while (std::getline(line_stream, value, ',')) {
            point.push_back(std::stod(value));
        }

        if (!points.empty() && point.size() != points[0].size()) {
            throw std::runtime_error("CSV file has rows with different sizes");
        }
        points.push_back(point);
    }
    return points;
}

void write_dendrogram_csv(const std::vector<MergeLine>& dendrogram,
                          const std::string& file_name) {
    std::ofstream output_file(file_name);
    if (!output_file.is_open()) {
        throw std::runtime_error("Could not write " + file_name);
    }

    output_file << "cl1,cl2,dist,new_size\n";
    for (const MergeLine& merge : dendrogram) {
        output_file << merge.first_cluster << ',' << merge.second_cluster << ','
                    << merge.distance << ',' << merge.merged_size << '\n';
    }
}

double euclidean_distance(const std::vector<double>& first_point,
                          const std::vector<double>& second_point) {
    double squared_sum = 0.0;
    for (size_t coordinate = 0; coordinate < first_point.size(); ++coordinate) {
        const double difference = first_point[coordinate] - second_point[coordinate];
        squared_sum += difference * difference;
    }
    return std::sqrt(squared_sum);
}

bool edge_is_better(const WeightedEdge& candidate,
                    const WeightedEdge& current_best) {
    if (candidate.distance != current_best.distance) {
        return candidate.distance < current_best.distance;
    }
    if (candidate.first_point != current_best.first_point) {
        return candidate.first_point < current_best.first_point;
    }
    return candidate.second_point < current_best.second_point;
}

void keep_cheapest_edge(std::vector<WeightedEdge>& cheapest_edges,
                        int component_id,
                        const WeightedEdge& edge) {
    if (cheapest_edges[component_id].first_point == -1 ||
        edge_is_better(edge, cheapest_edges[component_id])) {
        cheapest_edges[component_id] = edge;
    }
}

std::vector<WeightedEdge> build_mst_with_boruvka(
    const std::vector<std::vector<double>>& points,
    int number_of_threads) {
    const int number_of_points = static_cast<int>(points.size());
    if (number_of_points < 2) {
        throw std::runtime_error("Need at least two points");
    }
    if (number_of_threads < 1) {
        throw std::runtime_error("Need at least one thread");
    }

    number_of_threads = std::min(number_of_threads, number_of_points);

    Components components(number_of_points);
    int number_of_components = number_of_points;
    std::vector<WeightedEdge> mst_edges;
    mst_edges.reserve(number_of_points - 1);

    const WeightedEdge no_edge{-1, -1, std::numeric_limits<double>::infinity()};

    while (number_of_components > 1) {
        std::vector<std::vector<WeightedEdge>> thread_local_cheapest(
            number_of_threads, std::vector<WeightedEdge>(number_of_points, no_edge));

        std::vector<std::thread> workers;
        for (int thread_id = 0; thread_id < number_of_threads; ++thread_id) {
            const int first_i = (number_of_points * thread_id) / number_of_threads;
            const int after_last_i =
                (number_of_points * (thread_id + 1)) / number_of_threads;

            workers.emplace_back([&, thread_id, first_i, after_last_i]() {
                for (int first_point = first_i; first_point < after_last_i;
                     ++first_point) {
                    const int first_component =
                        components.find_without_compression(first_point);

                    for (int second_point = first_point + 1;
                         second_point < number_of_points;
                         ++second_point) {
                        const int second_component =
                            components.find_without_compression(second_point);
                        if (first_component == second_component) {
                            continue;
                        }

                        WeightedEdge edge{first_point, second_point,
                                          euclidean_distance(points[first_point],
                                                             points[second_point])};
                        keep_cheapest_edge(thread_local_cheapest[thread_id],
                                           first_component, edge);
                        keep_cheapest_edge(thread_local_cheapest[thread_id],
                                           second_component, edge);
                    }
                }
            });
        }

        for (std::thread& worker : workers) {
            worker.join();
        }

        std::vector<WeightedEdge> cheapest_edges(number_of_points, no_edge);
        for (const std::vector<WeightedEdge>& local_edges : thread_local_cheapest) {
            for (int component_id = 0; component_id < number_of_points;
                 ++component_id) {
                if (local_edges[component_id].first_point != -1) {
                    keep_cheapest_edge(cheapest_edges, component_id,
                                       local_edges[component_id]);
                }
            }
        }

        bool at_least_one_merge = false;
        for (int component_id = 0; component_id < number_of_points; ++component_id) {
            if (components.find_without_compression(component_id) != component_id ||
                cheapest_edges[component_id].first_point == -1) {
                continue;
            }

            const WeightedEdge edge = cheapest_edges[component_id];
            const int first_root = components.find(edge.first_point);
            const int second_root = components.find(edge.second_point);
            if (first_root == second_root) {
                continue;
            }

            components.merge_roots(first_root, second_root);
            mst_edges.push_back(edge);
            --number_of_components;
            at_least_one_merge = true;

            if (number_of_components == 1) {
                break;
            }
        }

        if (!at_least_one_merge) {
            throw std::runtime_error("MST construction got stuck");
        }
    }

    return mst_edges;
}

std::vector<MergeLine> mst_edges_to_dendrogram(
    int number_of_points,
    std::vector<WeightedEdge> mst_edges) {
    std::sort(mst_edges.begin(), mst_edges.end(),
              [](const WeightedEdge& first_edge, const WeightedEdge& second_edge) {
                  if (first_edge.distance != second_edge.distance) {
                      return first_edge.distance < second_edge.distance;
                  }
                  if (first_edge.first_point != second_edge.first_point) {
                      return first_edge.first_point < second_edge.first_point;
                  }
                  return first_edge.second_point < second_edge.second_point;
              });

    Components components(number_of_points);
    std::vector<int> cluster_name(number_of_points);
    std::iota(cluster_name.begin(), cluster_name.end(), 0);

    std::vector<MergeLine> dendrogram;
    dendrogram.reserve(number_of_points - 1);

    int next_cluster_name = number_of_points;
    for (const WeightedEdge& edge : mst_edges) {
        const int first_root = components.find(edge.first_point);
        const int second_root = components.find(edge.second_point);
        if (first_root == second_root) {
            continue;
        }

        const int merged_size =
            components.component_size[first_root] +
            components.component_size[second_root];
        dendrogram.push_back({cluster_name[first_root], cluster_name[second_root],
                              edge.distance, merged_size});

        const int new_root = components.merge_roots(first_root, second_root);
        cluster_name[new_root] = next_cluster_name++;
    }

    return dendrogram;
}

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "Usage: ./hac_single_mst <input.csv> <output.csv> <threads>\n";
        return 1;
    }

    try {
        const int number_of_threads = std::stoi(argv[3]);
        const std::vector<std::vector<double>> points = read_points_csv(argv[1]);

        const auto start_time = std::chrono::high_resolution_clock::now();
        std::vector<WeightedEdge> mst_edges =
            build_mst_with_boruvka(points, number_of_threads);
        const std::vector<MergeLine> dendrogram =
            mst_edges_to_dendrogram(static_cast<int>(points.size()), mst_edges);
        const auto end_time = std::chrono::high_resolution_clock::now();

        const double elapsed_ms =
            std::chrono::duration<double, std::milli>(end_time - start_time)
                .count();
        std::cout << "MST version finished in " << elapsed_ms << " ms using "
                  << number_of_threads << " thread(s)\n";
        write_dendrogram_csv(dendrogram, argv[2]);
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }

    return 0;
}
