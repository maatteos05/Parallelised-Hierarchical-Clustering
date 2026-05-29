import argparse
import csv
import math
import sys

TOL = 1e-5


def load_points(path: str) -> list[list[float]]:
    points = []
    with open(path, newline="") as csv_file:
        reader = csv.reader(csv_file)
        for row in reader:
            if row:
                points.append([float(value) for value in row])
    return points


def distance(first_point: list[float], second_point: list[float]) -> float:
    return math.sqrt(sum(
        (a - b) * (a - b)
        for a, b in zip(first_point, second_point)
    ))


def reference_average_link(points: list[list[float]]) -> list[tuple[float, int]]:
    number_of_points = len(points)
    maximum_cluster_count = 2 * number_of_points
    distances = [
        [0.0 for _ in range(maximum_cluster_count)]
        for _ in range(maximum_cluster_count)
    ]
    cluster_sizes = [0 for _ in range(maximum_cluster_count)]

    for first_point in range(number_of_points):
        cluster_sizes[first_point] = 1
        for second_point in range(first_point + 1, number_of_points):
            d = distance(points[first_point], points[second_point])
            distances[first_point][second_point] = d
            distances[second_point][first_point] = d

    active_clusters = list(range(number_of_points))
    reference_pairs = []
    next_cluster = number_of_points

    while len(active_clusters) > 1:
        best_distance = math.inf
        first_cluster = -1
        second_cluster = -1

        for i, cluster_i in enumerate(active_clusters):
            for cluster_j in active_clusters[i + 1:]:
                if distances[cluster_i][cluster_j] < best_distance:
                    best_distance = distances[cluster_i][cluster_j]
                    first_cluster = cluster_i
                    second_cluster = cluster_j

        new_cluster = next_cluster
        next_cluster += 1
        cluster_sizes[new_cluster] = (
            cluster_sizes[first_cluster] + cluster_sizes[second_cluster]
        )

        for other_cluster in active_clusters:
            if other_cluster in (first_cluster, second_cluster):
                continue

            weighted_distance = (
                cluster_sizes[first_cluster] * distances[first_cluster][other_cluster]
                + cluster_sizes[second_cluster] * distances[second_cluster][other_cluster]
            ) / cluster_sizes[new_cluster]
            distances[new_cluster][other_cluster] = weighted_distance
            distances[other_cluster][new_cluster] = weighted_distance

        active_clusters = [
            cluster
            for cluster in active_clusters
            if cluster not in (first_cluster, second_cluster)
        ]
        active_clusters.append(new_cluster)
        reference_pairs.append((best_distance, cluster_sizes[new_cluster]))

    return sorted(reference_pairs)


def load_cpp(path: str) -> list[tuple[float, int]]:
    pairs = []
    with open(path, newline="") as csv_file:
        reader = csv.DictReader(csv_file)
        for row in reader:
            pairs.append((float(row["dist"]), int(row["new_size"])))
    return sorted(pairs)


def compare(reference_pairs: list[tuple[float, int]],
            cpp_pairs: list[tuple[float, int]]) -> bool:
    if len(reference_pairs) != len(cpp_pairs):
        print(f"FAIL length mismatch: ref={len(reference_pairs)}, cpp={len(cpp_pairs)}")
        return False

    if all(abs(ref_dist - cpp_dist) <= TOL and ref_size == cpp_size
           for (ref_dist, ref_size), (cpp_dist, cpp_size)
           in zip(reference_pairs, cpp_pairs)):
        print(f"PASS {len(reference_pairs)} merges match (tol={TOL})")
        return True

    for index, (ref_row, cpp_row) in enumerate(zip(reference_pairs, cpp_pairs)):
        if abs(ref_row[0] - cpp_row[0]) > TOL or ref_row[1] != cpp_row[1]:
            print(f"FAIL first mismatch at sorted row {index}:")
            print(f"  ref: dist={ref_row[0]:.8f}, size={ref_row[1]}")
            print(f"  cpp: dist={cpp_row[0]:.8f}, size={cpp_row[1]}")
            break
    return False


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", required=True, help="Input CSV")
    parser.add_argument("--cpp_out", required=True, help="C++ dendrogram CSV")
    args = parser.parse_args()

    points = load_points(args.data)
    reference_pairs = reference_average_link(points)
    cpp_pairs = load_cpp(args.cpp_out)

    print(f"Points: {len(points)}")
    print(f"Merges: ref={len(reference_pairs)}, cpp={len(cpp_pairs)}")

    ok = compare(reference_pairs, cpp_pairs)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()