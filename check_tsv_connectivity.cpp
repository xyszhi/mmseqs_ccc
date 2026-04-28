#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
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
        "Check connectivity of a clusters TSV file (rep<TAB>member) using Union-Find.\n"
        "Reads the TSV directly without converting to METIS format first.\n"
        "\n"
        "Options:\n"
        "  --input   <path>  Input clusters TSV file (required)\n"
        "  -i        <path>  Alias for --input\n"
        "  --clique          Use full clique topology within each cluster\n"
        "                    (default: star topology, rep <-> each member)\n"
        "  --verbose         Print all component sizes (sorted descending, top 20)\n"
        "  --help            Show this help\n"
        "  -h                Show this help\n"
        "\n"
        "TSV format: each line is  <rep>TAB<member>\n"
        "  star   (default): rep <-> each member form edges\n"
        "  clique (--clique): all nodes within a cluster are pairwise connected\n"
        "\n"
        "Output:\n"
        "  Number of nodes, edges, connected components.\n"
        "  Size of the largest component.\n"
        "  With --verbose: full distribution of component sizes (top 20).\n"
        "\n"
        "Exit codes:\n"
        "  0 = connected (single component)\n"
        "  2 = disconnected (multiple components)\n"
        "  1 = error\n";
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
    bool full_clique = false;
    bool verbose = false;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--input" || arg == "-i") {
            if (i + 1 < argc) input_file = argv[++i];
        } else if (arg == "--clique") {
            full_clique = true;
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

    // ------------------------------------------------------------------
    // Pass 1: assign integer IDs to all node names
    // ------------------------------------------------------------------
    std::cerr << "[Pass 1] Scanning TSV and assigning node IDs...\n";

    std::unordered_map<std::string, int32_t> name_to_id;
    int32_t next_id = 1;

    {
        std::ifstream fin(input_file);
        if (!fin) {
            std::cerr << "Error: Cannot open file: " << input_file << "\n";
            return 1;
        }

        std::string line, rep, mem;
        while (std::getline(fin, line)) {
            if (line.empty() || line[0] == '%') continue;
            std::istringstream ss(line);
            if (!(ss >> rep >> mem)) continue;
            for (const std::string* s : {&rep, &mem}) {
                if (name_to_id.find(*s) == name_to_id.end()) {
                    name_to_id[*s] = next_id++;
                }
            }
        }
    }

    int32_t num_nodes = next_id - 1;
    std::cerr << "  " << num_nodes << " unique nodes found.\n";

    if (num_nodes == 0) {
        std::cerr << "Error: No nodes found in TSV file.\n";
        return 1;
    }

    // ------------------------------------------------------------------
    // Pass 2: stream TSV again, union edges
    // ------------------------------------------------------------------
    std::cerr << "[Pass 2] Building Union-Find from edges";
    if (full_clique) std::cerr << " (clique topology)";
    else             std::cerr << " (star topology)";
    std::cerr << "...\n";

    UnionFind uf(num_nodes);
    int64_t edges_processed = 0;

    {
        std::ifstream fin(input_file);
        if (!fin) {
            std::cerr << "Error: Cannot open file: " << input_file << "\n";
            return 1;
        }

        std::string line, rep, mem;

        if (full_clique) {
            // Clique mode: accumulate members per cluster, then union all pairs
            std::string cur_rep;
            std::vector<int32_t> cur_members; // includes rep_id

            auto flush_cluster = [&]() {
                if (cur_members.empty()) return;
                // union all pairs
                for (size_t i = 0; i < cur_members.size(); i++) {
                    for (size_t j = i + 1; j < cur_members.size(); j++) {
                        uf.unite(cur_members[i], cur_members[j]);
                        ++edges_processed;
                    }
                }
            };

            while (std::getline(fin, line)) {
                if (line.empty() || line[0] == '%') continue;
                std::istringstream ss(line);
                if (!(ss >> rep >> mem)) continue;

                if (rep != cur_rep) {
                    flush_cluster();
                    cur_rep = rep;
                    cur_members.clear();
                    cur_members.push_back(name_to_id.at(rep));
                }
                int32_t mem_id = name_to_id.at(mem);
                // avoid duplicates within cluster
                if (mem_id != cur_members[0]) {
                    cur_members.push_back(mem_id);
                }
            }
            flush_cluster();
        } else {
            // Star topology: one edge per line (rep <-> mem)
            while (std::getline(fin, line)) {
                if (line.empty() || line[0] == '%') continue;
                std::istringstream ss(line);
                if (!(ss >> rep >> mem)) continue;
                int32_t rep_id = name_to_id.at(rep);
                int32_t mem_id = name_to_id.at(mem);
                if (rep_id != mem_id) {
                    uf.unite(rep_id, mem_id);
                    ++edges_processed;
                }
            }
        }
    }

    std::cerr << "  " << edges_processed << " edges processed.\n";

    // ------------------------------------------------------------------
    // Count connected components
    // ------------------------------------------------------------------
    std::vector<int32_t> comp_size(num_nodes + 1, 0);
    for (int32_t i = 1; i <= num_nodes; i++) {
        comp_size[uf.find(i)]++;
    }

    int64_t num_components = 0;
    int64_t max_comp = 0;
    std::vector<int32_t> sizes;

    for (int32_t i = 1; i <= num_nodes; i++) {
        if (comp_size[i] > 0) {
            ++num_components;
            if (comp_size[i] > max_comp) max_comp = comp_size[i];
            if (verbose) sizes.push_back(comp_size[i]);
        }
    }

    // ------------------------------------------------------------------
    // Report
    // ------------------------------------------------------------------
    std::cout << "=== TSV Connectivity Report ===\n";
    std::cout << "Input file        : " << input_file << "\n";
    std::cout << "Topology          : " << (full_clique ? "clique" : "star") << "\n";
    std::cout << "Nodes             : " << num_nodes << "\n";
    std::cout << "Edges processed   : " << edges_processed << "\n";
    std::cout << "Connected components: " << num_components << "\n";
    std::cout << "Largest component : " << max_comp << " nodes ("
              << (100.0 * max_comp / num_nodes) << "%)\n";

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
            std::cout << "  #" << (i + 1) << ": " << sizes[i] << " nodes\n";
        }
        if ((int64_t)sizes.size() > 20) {
            std::cout << "  ... (" << (sizes.size() - 20) << " more components)\n";
        }
    }

    return (num_components == 1) ? 0 : 2;
}
