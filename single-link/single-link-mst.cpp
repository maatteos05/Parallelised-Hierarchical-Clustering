#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

struct merge_step {
    int left_cluster;
    int right_cluster;
    double distance;
    int new_cluster_size;
};

struct edge {
    int point_a;
    int point_b;
    double distance;
};

struct component_groups {
    std::vector<int> parent;
    std::vector<int> height;
    std::vector<int> group_size;

    component_groups(int nb_points) : parent(nb_points), height(nb_points, 0), group_size(nb_points, 1) {
        for (int point_id = 0; point_id < nb_points; point_id++) {
            parent[point_id] = point_id;
        }
    }

    int find_root(int point_id) {
        if (parent[point_id] != point_id) {
            parent[point_id] = find_root(parent[point_id]);
        }

        return parent[point_id];
    }

    int find_root_read_only(int point_id) const {
        while (parent[point_id] != point_id) {
            point_id = parent[point_id];
        }

        return point_id;
    }

    int attach_two_roots(int root_a, int root_b) {
        if (height[root_a] < height[root_b]) {
            int temporary_root = root_a;
            root_a = root_b;
            root_b = temporary_root;
        }

        parent[root_b] = root_a;
        group_size[root_a] = group_size[root_a] + group_size[root_b];

        if (height[root_a] == height[root_b]) {
            height[root_a]++;
        }

        return root_a;
    }
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
    std::vector<std::vector<double>> points;
    std::ifstream input_file(input_path);

    if (!input_file.is_open()) {
        throw std::runtime_error("Could not open " + input_path);
    }

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

bool edge_is_better(const edge& candidate_edge, const edge& current_edge) {
    if (candidate_edge.distance < current_edge.distance) {
        return true;
    }

    if (candidate_edge.distance > current_edge.distance) {
        return false;
    }

    if (candidate_edge.point_a < current_edge.point_a) {
        return true;
    }

    if (candidate_edge.point_a > current_edge.point_a) {
        return false;
    }

    return candidate_edge.point_b < current_edge.point_b;
}

bool edge_goes_first(const edge& edge_a, const edge& edge_b) {
    return edge_is_better(edge_a, edge_b);
}

void sort_edges_by_distance(std::vector<edge>& edges) {
    for (int i = 0; i < static_cast<int>(edges.size()); i++) {
        int best_position = i;

        for (int j = i + 1; j < static_cast<int>(edges.size()); j++) {
            if (edge_goes_first(edges[j], edges[best_position])) {
                best_position = j;
            }
        }

        if (best_position != i) {
            edge temporary_edge = edges[i];
            edges[i] = edges[best_position];
            edges[best_position] = temporary_edge;
        }
    }
}

void keep_if_better(std::vector<edge>& best_edges, int group_id, const edge& candidate_edge) {
    if (best_edges[group_id].point_a == -1) {
        best_edges[group_id] = candidate_edge;
    } else if (edge_is_better(candidate_edge, best_edges[group_id])) {
        best_edges[group_id] = candidate_edge;
    }
}

edge make_empty_edge() {
    edge empty_edge;
    empty_edge.point_a = -1;
    empty_edge.point_b = -1;
    empty_edge.distance = std::numeric_limits<double>::infinity();
    return empty_edge;
}

void scan_points(const std::vector<std::vector<double>>* points, const component_groups* groups, int start_point, int end_point, int thread_id, std::vector<std::vector<edge>>* local_best) {
    int nb_points = static_cast<int>(points->size());

    for (int point_a = start_point; point_a < end_point; point_a++) {
        int group_a = groups->find_root_read_only(point_a);

        for (int point_b = point_a + 1; point_b < nb_points; point_b++) {
            int group_b = groups->find_root_read_only(point_b);

            if (group_a != group_b) {
                edge current_edge;
                current_edge.point_a = point_a;
                current_edge.point_b = point_b;
                current_edge.distance = distance_between_points((*points)[point_a], (*points)[point_b]);

                keep_if_better((*local_best)[thread_id], group_a, current_edge);
                keep_if_better((*local_best)[thread_id], group_b, current_edge);
            }
        }
    }
}

std::vector<edge> build_mst_with_boruvka(const std::vector<std::vector<double>>& points, int nb_threads) {
    int nb_points = static_cast<int>(points.size());

    if (nb_points < 2) {
        throw std::runtime_error("Need at least two points");
    }

    if (nb_threads < 1) {
        throw std::runtime_error("Need at least one thread");
    }

    if (nb_threads > nb_points) {
        nb_threads = nb_points;
    }

    edge empty_edge = make_empty_edge();
    component_groups groups(nb_points);
    int nb_groups_left = nb_points;
    std::vector<edge> mst_edges;

    while (nb_groups_left > 1) {
        std::vector<std::vector<edge>> local_best(nb_threads, std::vector<edge>(nb_points, empty_edge));
        std::vector<std::thread> workers;

        for (int thread_id = 0; thread_id < nb_threads; thread_id++) {
            int start_point = (nb_points * thread_id) / nb_threads;
            int end_point = (nb_points * (thread_id + 1)) / nb_threads;

            workers.push_back(std::thread(scan_points, &points, &groups, start_point, end_point, thread_id, &local_best));
        }

        for (int thread_id = 0; thread_id < static_cast<int>(workers.size()); thread_id++) {
            workers[thread_id].join();
        }

        std::vector<edge> best_edge(nb_points, empty_edge);

        for (int thread_id = 0; thread_id < nb_threads; thread_id++) {
            for (int group_id = 0; group_id < nb_points; group_id++) {
                if (local_best[thread_id][group_id].point_a != -1) {
                    keep_if_better(best_edge, group_id, local_best[thread_id][group_id]);
                }
            }
        }

        bool merged_something = false;

        for (int group_id = 0; group_id < nb_points; group_id++) {
            if (groups.find_root_read_only(group_id) != group_id) {
                continue;
            }

            if (best_edge[group_id].point_a != -1) {
                edge chosen_edge = best_edge[group_id];
                int root_a = groups.find_root(chosen_edge.point_a);
                int root_b = groups.find_root(chosen_edge.point_b);

                if (root_a != root_b) {
                    groups.attach_two_roots(root_a, root_b);
                    mst_edges.push_back(chosen_edge);
                    nb_groups_left--;
                    merged_something = true;
                }
            }

            if (nb_groups_left == 1) {
                break;
            }
        }

        if (!merged_something) {
            throw std::runtime_error("MST construction got stuck");
        }
    }

    return mst_edges;
}

std::vector<merge_step> mst_edges_to_dendrogram(int nb_points, std::vector<edge> mst_edges) {
    sort_edges_by_distance(mst_edges);

    component_groups groups(nb_points);
    std::vector<int> cluster_name(nb_points);

    for (int point_id = 0; point_id < nb_points; point_id++) {
        cluster_name[point_id] = point_id;
    }

    std::vector<merge_step> dendrogram;
    int next_cluster_name = nb_points;

    for (int i = 0; i < static_cast<int>(mst_edges.size()); i++) {
        edge current_edge = mst_edges[i];
        int root_a = groups.find_root(current_edge.point_a);
        int root_b = groups.find_root(current_edge.point_b);

        if (root_a != root_b) {
            merge_step new_line;
            new_line.left_cluster = cluster_name[root_a];
            new_line.right_cluster = cluster_name[root_b];
            new_line.distance = current_edge.distance;
            new_line.new_cluster_size = groups.group_size[root_a] + groups.group_size[root_b];
            dendrogram.push_back(new_line);

            int new_root = groups.attach_two_roots(root_a, root_b);
            cluster_name[new_root] = next_cluster_name;
            next_cluster_name++;
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
    if (argc != 4) {
        std::cerr << "Usage: ./single_link_mst <input.csv> <output.csv> <threads>\n";
        return 1;
    }

    try {
        int nb_threads = std::stoi(argv[3]);
        std::vector<std::vector<double>> points = read_points_csv(argv[1]);

        std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();
        std::vector<edge> mst_edges = build_mst_with_boruvka(points, nb_threads);
        std::vector<merge_step> dendrogram = mst_edges_to_dendrogram(static_cast<int>(points.size()), mst_edges);
        std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();

        double elapsed_ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();
        std::cout << "MST version finished in " << elapsed_ms << " ms using " << nb_threads << " thread(s)\n";
        write_dendrogram_csv(dendrogram, argv[2]);
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }

    return 0;
}
