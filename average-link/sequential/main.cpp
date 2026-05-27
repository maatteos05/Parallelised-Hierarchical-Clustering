#include "average-link-seq.hpp"
#include <chrono>
#include <iostream>
#include <stdexcept>

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: hac_seq <input.csv> <output.csv>\n";
        return 1;
    }

    const std::string input_path  = argv[1];
    const std::string output_path = argv[2];

    try {
        // ── Load data ─────────────────────────────────────────────────────
        auto data = load_csv(input_path);
        std::cout << "Loaded " << data.size() << " points"
                  << " (dim=" << (data.empty() ? 0 : data[0].size()) << ")\n";

        // ── Run HAC ───────────────────────────────────────────────────────
        auto t0 = std::chrono::high_resolution_clock::now();

        auto dendrogram = hac_average_link(data);

        auto t1  = std::chrono::high_resolution_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

        std::cout << "HAC done in " << ms << " ms  ("
                  << dendrogram.size() << " merges)\n";

        // ── Save result ───────────────────────────────────────────────────
        save_dendrogram(dendrogram, output_path);
        std::cout << "Dendrogram saved to " << output_path << "\n";

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}