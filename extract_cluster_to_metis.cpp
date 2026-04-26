#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <set>
#include <algorithm>

// ---------------------------------------------------------------------------
// Usage
// ---------------------------------------------------------------------------
static void print_usage(const char* prog) {
    std::cerr <<
        "Usage: " << prog << " [OPTIONS]\n"
        "\n"
        "Convert a clusters TSV file (rep<TAB>member) to METIS graph format.\n"
        "\n"
        "Options:\n"
        "  --input  <path>  Input clusters TSV file (default: test/filtered_db_clusters.tsv)\n"
        "  -i       <path>  Alias for --input\n"
        "  --output <path>  Output METIS .graph file (default: <input>.graph)\n"
        "  -o       <path>  Alias for --output\n"
        "  --clique         Build full clique within each cluster (default: star topology)\n"
        "  --help           Show this help\n"
        "  -h               Show this help\n"
        "\n"
        "Topology modes:\n"
        "  star   (default): representative <-> each member (fewer edges)\n"
        "  clique (--clique): all pairs within cluster connected (more edges)\n"
        "\n"
        "Output format: METIS unweighted undirected graph (.graph)\n"
        "  Line 1: <num_vertices> <num_edges>\n"
        "  Line i: space-separated neighbor IDs (1-based) of vertex i\n";
}

int main(int argc, char* argv[]) {
    std::string input_file = "test/filtered_db_clusters.tsv";
    std::string output_file;
    bool full_clique = false;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--input" || arg == "-i") {
            if (i + 1 < argc) input_file = argv[++i];
        } else if (arg == "--output" || arg == "-o") {
            if (i + 1 < argc) output_file = argv[++i];
        } else if (arg == "--clique") {
            full_clique = true;
        } else if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            return 0;
        } else {
            std::cerr << "Unknown argument: " << arg << "\n";
            print_usage(argv[0]);
            return 1;
        }
    }

    if (output_file.empty()) {
        output_file = input_file + ".graph";
    }

    // Step 1: read TSV, assign integer IDs, collect clusters
    std::unordered_map<std::string, int> name_to_id;
    std::unordered_map<std::string, std::vector<std::string>> clusters; // rep -> members

    {
        std::ifstream fin(input_file);
        if (!fin) {
            std::cerr << "Error: cannot open input file: " << input_file << "\n";
            return 1;
        }
        std::string line, rep, mem;
        while (std::getline(fin, line)) {
            if (line.empty() || line[0] == '%') continue;
            std::istringstream ss(line);
            if (!(ss >> rep >> mem)) continue;
            for (const auto& s : {rep, mem}) {
                if (name_to_id.find(s) == name_to_id.end()) {
                    int id = (int)name_to_id.size() + 1; // 1-based
                    name_to_id[s] = id;
                }
            }
            clusters[rep].push_back(mem);
        }
    }

    int num_vertices = (int)name_to_id.size();

    // Step 2: build adjacency lists (undirected, no self-loops, no duplicates)
    std::vector<std::set<int>> adj(num_vertices + 1);

    for (auto& kv : clusters) {
        const std::string& rep = kv.first;
        const std::vector<std::string>& members = kv.second;
        int rep_id = name_to_id[rep];

        if (full_clique) {
            // full clique: all pairs within cluster
            for (size_t i = 0; i < members.size(); i++) {
                int u = name_to_id[members[i]];
                if (u != rep_id) {
                    adj[rep_id].insert(u);
                    adj[u].insert(rep_id);
                }
                for (size_t j = i + 1; j < members.size(); j++) {
                    int v = name_to_id[members[j]];
                    if (u != v) {
                        adj[u].insert(v);
                        adj[v].insert(u);
                    }
                }
            }
        } else {
            // star: representative <-> each member
            for (const auto& mem : members) {
                int mem_id = name_to_id[mem];
                if (mem_id == rep_id) continue; // skip self-loop
                adj[rep_id].insert(mem_id);
                adj[mem_id].insert(rep_id);
            }
        }
    }

    // count edges (each undirected edge counted once)
    long long num_edges = 0;
    for (int i = 1; i <= num_vertices; i++) {
        num_edges += (long long)adj[i].size();
    }
    num_edges /= 2;

    // Step 3: write METIS format
    std::ofstream fout(output_file);
    if (!fout) {
        std::cerr << "Error: cannot write output file: " << output_file << "\n";
        return 1;
    }

    fout << num_vertices << " " << num_edges << "\n";
    for (int i = 1; i <= num_vertices; i++) {
        bool first = true;
        for (int nb : adj[i]) {
            if (!first) fout << " ";
            fout << nb;
            first = false;
        }
        fout << "\n";
    }

    std::cerr << "Done: " << num_vertices << " vertices, " << num_edges << " edges"
              << " (" << (full_clique ? "clique" : "star") << " topology)"
              << " -> " << output_file << "\n";
    return 0;
}
