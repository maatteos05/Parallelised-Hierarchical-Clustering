#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

struct merge_step {
    int left_cluster;
    int right_cluster;
    double distance;
    int new_cluster_size;
};

double distance_between_points(const std::vector<double>& point_a, const std::vector<double>& point_b) {
    double sum_of_squares = 0.0;

    for (int coord = 0; coord < static_cast<int>(point_a.size()); coord++) {
        double difference = point_a[coord] - point_b[coord];
        sum_of_squares = sum_of_squares + difference * difference;
    }

    return std::sqrt(sum_of_squares);
}

std::vector<std::vector<double>> read_points_csv(const std::string& input_path) {
    std::ifstream input_file(input_path);
    if (!input_file.is_open()) {
        throw std::runtime_error("Could not open " + input_path);
    }

    std::vector<std::vector<double>> points;
    std::string line;

    while (std::getline(input_file, line)) {
        if (line.size() == 0) {
            continue;
        }

        std::vector<double> current_point;
        std::stringstream line_stream(line);
        std::string cell;

        while (std::getline(line_stream, cell, ',')) {
            current_point.push_back(std::stod(cell));
        }

        if (points.size() > 0 && current_point.size() != points[0].size()) {
            throw std::runtime_error("CSV file has rows with different sizes");
        }

        points.push_back(current_point);
    }

    return points;
}

int closest_alive_cluster(int cluster_to_check, const std::vector<int>& alive, const std::vector<std::vector<double>>& distances) {
    int closest_cluster = -1;
    double closest_distance = std::numeric_limits<double>::infinity();

    for (int index = 0; index < static_cast<int>(alive.size()); index++) {
        int other_cluster = alive[index];

        if (other_cluster != cluster_to_check) {
            double current_distance = distances[cluster_to_check][other_cluster];

            if (current_distance < closest_distance) {
                closest_distance = current_distance;
                closest_cluster = other_cluster;
            }
        }
    }

    return closest_cluster;
}

std::vector<merge_step> single_link_baseline(const std::vector<std::vector<double>>& points) {
    int nb_points = static_cast<int>(points.size());

    if (nb_points < 2) {
        throw std::runtime_error("Need at least two points");
    }

    int max_nb_clusters = 2 * nb_points;
    std::vector<int> cluster_size(max_nb_clusters, 0);
    std::vector<int> alive;

    for (int point_id = 0; point_id < nb_points; point_id++) {
        alive.push_back(point_id);
        cluster_size[point_id] = 1;
    }

    std::vector<std::vector<double>> distances(max_nb_clusters, std::vector<double>(max_nb_clusters, 0.0));

    for (int point_a = 0; point_a < nb_points; point_a++) {
        for (int point_b = point_a + 1; point_b < nb_points; point_b++) {
            double distance = distance_between_points(points[point_a], points[point_b]);
            distances[point_a][point_b] = distance;
            distances[point_b][point_a] = distance;
        }
    }

    std::vector<int> closest_cluster(max_nb_clusters, -1);
    for (int i = 0; i < static_cast<int>(alive.size()); i++) {
        int cluster_id = alive[i];
        closest_cluster[cluster_id] = closest_alive_cluster(cluster_id, alive, distances);
    }

    std::vector<merge_step> dendrogram;
    int next_cluster_id = nb_points;

    while (alive.size() > 1) {
        int cluster_a = -1;
        int cluster_b = -1;
        double best_distance = std::numeric_limits<double>::infinity();

        for (int i = 0; i < static_cast<int>(alive.size()); i++) {
            int current_cluster = alive[i];
            int nearest = closest_cluster[current_cluster];

            if (nearest != -1) {
                double current_distance = distances[current_cluster][nearest];

                if (current_distance < best_distance) {
                    best_distance = current_distance;
                    cluster_a = current_cluster;
                    cluster_b = nearest;
                }
            }
        }

        int new_cluster = next_cluster_id;
        next_cluster_id++;
        cluster_size[new_cluster] = cluster_size[cluster_a] + cluster_size[cluster_b];

        merge_step new_line;
        new_line.left_cluster = cluster_a;
        new_line.right_cluster = cluster_b;
        new_line.distance = best_distance;
        new_line.new_cluster_size = cluster_size[new_cluster];
        dendrogram.push_back(new_line);

        std::vector<int> next_alive_clusters;
        for (int i = 0; i < static_cast<int>(alive.size()); i++) {
            int other_cluster = alive[i];

            if (other_cluster == cluster_a || other_cluster == cluster_b) {
                continue;
            }

            double distance_from_a = distances[cluster_a][other_cluster];
            double distance_from_b = distances[cluster_b][other_cluster];
            double new_distance = distance_from_a;

            if (distance_from_b < new_distance) {
                new_distance = distance_from_b;
            }

            distances[new_cluster][other_cluster] = new_distance;
            distances[other_cluster][new_cluster] = new_distance;
            next_alive_clusters.push_back(other_cluster);
        }

        next_alive_clusters.push_back(new_cluster);
        alive = next_alive_clusters;

        closest_cluster[new_cluster] = closest_alive_cluster(new_cluster, alive, distances);

        for (int i = 0; i < static_cast<int>(alive.size()); i++) {
            int current_cluster = alive[i];

            if (current_cluster != new_cluster) {
                int old_closest = closest_cluster[current_cluster];

                if (old_closest == cluster_a || old_closest == cluster_b) {
                    closest_cluster[current_cluster] = closest_alive_cluster(current_cluster, alive, distances);
                } else {
                    double old_distance = distances[current_cluster][old_closest];
                    double distance_to_new = distances[current_cluster][new_cluster];

                    if (distance_to_new < old_distance) {
                        closest_cluster[current_cluster] = new_cluster;
                    }
                }
            }
        }
    }

    return dendrogram;
}

void write_dendrogram_csv(const std::vector<merge_step>& dendrogram, const std::string& output_path) {
    std::ofstream output_file(output_path);
    if (!output_file.is_open()) {
        throw std::runtime_error("Could not write " + output_path);
    }

    output_file << "cl1,cl2,dist,new_size\n";
    for (int i = 0; i < static_cast<int>(dendrogram.size()); i++) {
        output_file << dendrogram[i].left_cluster << ',' << dendrogram[i].right_cluster << ',' << dendrogram[i].distance << ',' << dendrogram[i].new_cluster_size << '\n';
    }
}

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: ./hac_single_baseline <input.csv> <output.csv>\n";
        return 1;
    }

    try {
        std::vector<std::vector<double>> points = read_points_csv(argv[1]);

        std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();
        std::vector<merge_step> dendrogram = single_link_baseline(points);
        std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();

        double elapsed_ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();

        std::cout << "Baseline finished in " << elapsed_ms << " ms\n";
        write_dendrogram_csv(dendrogram, argv[2]);
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }

    return 0;
}
