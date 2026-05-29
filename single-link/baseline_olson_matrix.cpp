#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

struct MergeLine {
    int first_cluster;
    int second_cluster;
    double distance;
    int merged_size;
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

int nearest_active_cluster(
    int cluster_id,
    const std::vector<int>& active_clusters,
    const std::vector<std::vector<double>>& distance_matrix) {
    int nearest_cluster = -1;
    double nearest_distance = std::numeric_limits<double>::infinity();

    for (int other_cluster : active_clusters) {
        if (other_cluster == cluster_id) {
            continue;
        }

        const double current_distance = distance_matrix[cluster_id][other_cluster];
        if (current_distance < nearest_distance) {
            nearest_distance = current_distance;
            nearest_cluster = other_cluster;
        }
    }

    return nearest_cluster;
}

std::vector<MergeLine> single_link_baseline_olson_matrix(
    const std::vector<std::vector<double>>& points) {
    const int number_of_points = static_cast<int>(points.size());
    if (number_of_points < 2) {
        throw std::runtime_error("Need at least two points");
    }

    const int maximum_cluster_count = 2 * number_of_points;
    std::vector<std::vector<double>> distance_matrix(
        maximum_cluster_count,
        std::vector<double>(maximum_cluster_count, 0.0));

    for (int first_point = 0; first_point < number_of_points; ++first_point) {
        for (int second_point = first_point + 1; second_point < number_of_points;
             ++second_point) {
            const double distance =
                euclidean_distance(points[first_point], points[second_point]);
            distance_matrix[first_point][second_point] = distance;
            distance_matrix[second_point][first_point] = distance;
        }
    }

    std::vector<int> active_clusters;
    std::vector<int> cluster_sizes(maximum_cluster_count, 0);
    for (int point_id = 0; point_id < number_of_points; ++point_id) {
        active_clusters.push_back(point_id);
        cluster_sizes[point_id] = 1;
    }

    std::vector<int> nearest_neighbor(maximum_cluster_count, -1);
    for (int cluster_id : active_clusters) {
        nearest_neighbor[cluster_id] =
            nearest_active_cluster(cluster_id, active_clusters, distance_matrix);
    }

    std::vector<MergeLine> dendrogram;
    dendrogram.reserve(number_of_points - 1);

    int next_cluster_id = number_of_points;
    while (active_clusters.size() > 1) {
        int first_cluster_to_merge = -1;
        int second_cluster_to_merge = -1;
        double best_distance = std::numeric_limits<double>::infinity();

        for (int cluster_id : active_clusters) {
            const int candidate_neighbor = nearest_neighbor[cluster_id];
            if (candidate_neighbor == -1) {
                continue;
            }

            const double candidate_distance =
                distance_matrix[cluster_id][candidate_neighbor];
            if (candidate_distance < best_distance) {
                best_distance = candidate_distance;
                first_cluster_to_merge = cluster_id;
                second_cluster_to_merge = candidate_neighbor;
            }
        }

        const int new_cluster_id = next_cluster_id++;
        cluster_sizes[new_cluster_id] =
            cluster_sizes[first_cluster_to_merge] +
            cluster_sizes[second_cluster_to_merge];

        std::vector<int> new_active_clusters;
        new_active_clusters.reserve(active_clusters.size() - 1);
        for (int cluster_id : active_clusters) {
            if (cluster_id == first_cluster_to_merge ||
                cluster_id == second_cluster_to_merge) {
                continue;
            }

            const double new_distance = std::min(
                distance_matrix[first_cluster_to_merge][cluster_id],
                distance_matrix[second_cluster_to_merge][cluster_id]);
            distance_matrix[new_cluster_id][cluster_id] = new_distance;
            distance_matrix[cluster_id][new_cluster_id] = new_distance;
            new_active_clusters.push_back(cluster_id);
        }
        new_active_clusters.push_back(new_cluster_id);
        active_clusters = new_active_clusters;

        dendrogram.push_back({first_cluster_to_merge, second_cluster_to_merge,
                              best_distance, cluster_sizes[new_cluster_id]});

        nearest_neighbor[new_cluster_id] =
            nearest_active_cluster(new_cluster_id, active_clusters, distance_matrix);

        for (int cluster_id : active_clusters) {
            if (cluster_id == new_cluster_id) {
                continue;
            }

            const int old_neighbor = nearest_neighbor[cluster_id];
            if (old_neighbor == first_cluster_to_merge ||
                old_neighbor == second_cluster_to_merge) {
                nearest_neighbor[cluster_id] =
                    nearest_active_cluster(cluster_id, active_clusters,
                                           distance_matrix);
            } else if (distance_matrix[cluster_id][new_cluster_id] <
                       distance_matrix[cluster_id][old_neighbor]) {
                nearest_neighbor[cluster_id] = new_cluster_id;
            }
        }
    }

    return dendrogram;
}

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: ./hac_single_baseline <input.csv> <output.csv>\n";
        return 1;
    }

    try {
        const std::vector<std::vector<double>> points = read_points_csv(argv[1]);
        const auto start_time = std::chrono::high_resolution_clock::now();
        const std::vector<MergeLine> dendrogram =
            single_link_baseline_olson_matrix(points);
        const auto end_time = std::chrono::high_resolution_clock::now();

        const double elapsed_ms =
            std::chrono::duration<double, std::milli>(end_time - start_time)
                .count();
        std::cout << "Baseline finished in " << elapsed_ms << " ms\n";
        write_dendrogram_csv(dendrogram, argv[2]);
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }

    return 0;
}
