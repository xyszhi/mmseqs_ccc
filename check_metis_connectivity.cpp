#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <stdexcept>

// ---------------------------------------------------------------------------
// Usage
// ---------------------------------------------------------------------------
static void print_usage(const char* prog) {
    std::cerr <<
        "Usage: " << prog << " [OPTIONS]\n"
        "\n"
        "Check connectivity of a METIS graph file using Union-Find (并查集).\n"
        "Reads the METIS .graph file and reports connected components.\n"
        "\n"
        "Options:\n"
        "  --input  <path>  Input METIS .graph file (required)\n"
        "  -i       <path>  Alias for --input\n"
        "  --verbose        Print all component sizes (sorted descending)\n"
        "  --help           Show this help\n"
        "  -h               Show this help\n"
        "\n"
        "Output:\n"
        "  Number of vertices, edges, connected components.\n"
        "  Size of the largest component.\n"
        "  With --verbose: full distribution of component sizes.\n";
}

// ---------------------------------------------------------------------------
// Union-Find (path compression + union by rank)
// ---------------------------------------------------------------------------
struct UnionFind {
    std::vector<int32_t> parent;
    std::vector<int32_t> rank_;

    explicit UnionFind(int32_t n) : parent(n + 1), rank_(n + 1, 0) {
        for (int32_t i = 0; i <= n; i++) parent[i] = i;
    }

    int32_t find(int32_t x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]]; // path halving
            x = parent[x];
        }
        return x;
    }

    // Returns true if they were in different sets (i.e., a merge happened)
    bool unite(int32_t a, int32_t b) {
        a = find(a); b = find(b);
        if (a == b) return false;
        if (rank_[a] < rank_[b]) std::swap(a, b);
        parent[b] = a;
        if (rank_[a] == rank_[b]) rank_[a]++;
        return true;
    }
};

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    std::string input_file;
    bool verbose = false;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--input" || arg == "-i") {
            if (i + 1 < argc) input_file = argv[++i];
        } else if (arg == "--verbose") {
            verbose = true;
        } else if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            return 0;
        } else {
            std::cerr << "Unknown argument: " << arg << "\n";
            print_usage(argv[0]);
            return 1;
        }
    }

    if (input_file.empty()) {
        std::cerr << "Error: --input is required.\n\n";
        print_usage(argv[0]);
        return 1;
    }

    std::ifstream fin(input_file);
    if (!fin) {
        std::cerr << "Error: Cannot open file: " << input_file << "\n";
        return 1;
    }

    // ------------------------------------------------------------------
    // Read METIS header: first non-comment line is "<num_vertices> <num_edges> [fmt] ..."
    // ------------------------------------------------------------------
    std::string line;
    int64_t num_vertices = 0, num_edges_header = 0;

    while (std::getline(fin, line)) {
        if (line.empty() || line[0] == '%') continue;
        std::istringstream ss(line);
        ss >> num_vertices >> num_edges_header;
        break;
    }

    if (num_vertices <= 0) {
        std::cerr << "Error: Could not read valid METIS header.\n";
        return 1;
    }

    std::cerr << "Graph: " << num_vertices << " vertices, "
              << num_edges_header << " edges (from header)\n";
    std::cerr << "Building Union-Find structure...\n";

    UnionFind uf((int32_t)num_vertices);

    // ------------------------------------------------------------------
    // Read adjacency lines (one per vertex, 1-based neighbor IDs)
    // ------------------------------------------------------------------
    int64_t vertex = 0;
    int64_t edges_seen = 0;

    while (std::getline(fin, line)) {
        if (line.empty() || line[0] == '%') continue;
        ++vertex;
        if (vertex > num_vertices) break;

        std::istringstream ss(line);
        int32_t nb;
        while (ss >> nb) {
            if (nb < 1 || nb > (int32_t)num_vertices) continue;
            // Only process each undirected edge once (u < v)
            if ((int32_t)vertex < nb) {
                uf.unite((int32_t)vertex, nb);
                ++edges_seen;
            }
        }
    }

    std::cerr << "Edges processed (u<v): " << edges_seen << "\n";

    // ------------------------------------------------------------------
    // Count connected components
    // ------------------------------------------------------------------
    // component size: root -> count
    std::vector<int32_t> comp_size(num_vertices + 1, 0);
    for (int32_t i = 1; i <= (int32_t)num_vertices; i++) {
        comp_size[uf.find(i)]++;
    }

    int64_t num_components = 0;
    int64_t max_comp = 0;
    std::vector<int32_t> sizes;

    for (int32_t i = 1; i <= (int32_t)num_vertices; i++) {
        if (comp_size[i] > 0) {
            ++num_components;
            if (comp_size[i] > max_comp) max_comp = comp_size[i];
            if (verbose) sizes.push_back(comp_size[i]);
        }
    }

    // ------------------------------------------------------------------
    // Report
    // ------------------------------------------------------------------
    std::cout << "=== Connectivity Report ===\n";
    std::cout << "Vertices          : " << num_vertices << "\n";
    std::cout << "Edges (header)    : " << num_edges_header << "\n";
    std::cout << "Edges (processed) : " << edges_seen << "\n";
    std::cout << "Connected components: " << num_components << "\n";
    std::cout << "Largest component : " << max_comp << " vertices ("
              << (100.0 * max_comp / num_vertices) << "%)\n";

    if (num_components == 1) {
        std::cout << "Result: CONNECTED (single component)\n";
    } else {
        std::cout << "Result: DISCONNECTED (" << num_components << " components)\n";
    }

    if (verbose && !sizes.empty()) {
        std::sort(sizes.begin(), sizes.end(), std::greater<int32_t>());
        std::cout << "\n--- Component size distribution (top 20) ---\n";
        int show = (int)std::min((int64_t)20, (int64_t)sizes.size());
        for (int i = 0; i < show; i++) {
            std::cout << "  #" << (i + 1) << ": " << sizes[i] << " vertices\n";
        }
        if ((int64_t)sizes.size() > 20) {
            std::cout << "  ... (" << (sizes.size() - 20) << " more components)\n";
        }
    }

    return (num_components == 1) ? 0 : 2; // exit 0=connected, 2=disconnected
}
