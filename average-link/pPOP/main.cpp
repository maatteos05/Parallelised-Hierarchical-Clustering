#include "average-link-pPOP.hpp"
#include <chrono>
#include <iostream>
#include <string>
#include <cmath>

int main(int argc, char *argv[]) {
  if (argc < 4) {
    std::cerr << "Usage: hac_ppop <input.csv> <output.csv>"
                 " <n_threads> [initial_n_cells_per_dim]\n";
    return 1;
  }

  std::string input_path = argv[1];
  std::string output_path = argv[2];

  int n_threads = std::stoi(argv[3]);

  try {
    // ── Load data ─────────────────────────────────────────────────────
    auto data = load_csv(input_path);
    int N = (int)data.size();

    int n_cells_per_dim = std::max(1, (int)std::sqrt(N / 10.0));
    if (argc >= 5) {
      n_cells_per_dim = std::stoi(argv[4]);
    }

    std::cout << "Loaded " << data.size() << " points"
              << " (dim=" << (data.empty() ? 0 : data[0].size()) << ")\n";
    std::cout << "Running pPOP average-link with parameters:\n"
              << "  threads       = " << n_threads << "\n"
              << "  initial_cells = " << n_cells_per_dim << "\n";

    // ── Run HAC pPOP ──────────────────────────────────────────────────
    auto t0 = std::chrono::high_resolution_clock::now();

    auto dendrogram = hac_average_link_ppop(data, n_threads, n_cells_per_dim);

    auto t1 = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    std::cout << "HAC done in " << ms << " ms  (" << dendrogram.size()
              << " merges)\n";

    // ── Save result ───────────────────────────────────────────────────
    save_dendrogram(dendrogram, output_path);
    std::cout << "Dendrogram saved to " << output_path << "\n";

  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 1;
  }

  return 0;
}
