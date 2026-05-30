#include "average-link-naive.hpp"
#include <chrono>
#include <iostream>
#include <string>

int main(int argc, char *argv[]) {
  if (argc < 3) {
    std::cerr << "Usage: hac_naive <input.csv> <output.csv> [n_threads]\n";
    std::cerr << "Defaults: n_threads=4\n";
    return 1;
  }

  const std::string input_path = argv[1];
  const std::string output_path = argv[2];

  int n_threads = 4;
  if (argc >= 4)
    n_threads = std::stoi(argv[3]);

  try {
    auto data = load_csv(input_path);
    std::cout << "Loaded " << data.size() << " points"
              << " (dim=" << (data.empty() ? 0 : data[0].size()) << ")\n";
    std::cout << "Running naive parallel average-link with parameters:\n"
              << "  n_threads = " << n_threads << "\n";

    auto t0 = std::chrono::high_resolution_clock::now();

    auto dendrogram = hac_average_link_naive(data, n_threads);

    auto t1 = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    std::cout << "HAC done in " << ms << " ms  (" << dendrogram.size()
              << " merges)\n";

    save_dendrogram(dendrogram, output_path);
    std::cout << "Dendrogram saved to " << output_path << "\n";

  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 1;
  }

  return 0;
}
