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
    std::string map_file = output_file + ".id_map";

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

    int num_vertices_raw = (int)name_to_id.size();

    // Step 2: build adjacency lists (undirected, no self-loops, no duplicates)
    std::vector<std::set<int>> adj(num_vertices_raw + 1);

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

    // Step 3: compact renumbering — skip isolated vertices (no edges)
    // old_id (1-based) -> new_id (1-based), 0 means isolated/dropped
    std::vector<int> new_id(num_vertices_raw + 1, 0);
    int num_vertices = 0;
    for (int i = 1; i <= num_vertices_raw; i++) {
        if (!adj[i].empty()) {
            new_id[i] = ++num_vertices;
        }
    }

    // remap adjacency lists to new IDs
    std::vector<std::vector<int>> new_adj(num_vertices + 1);
    for (int i = 1; i <= num_vertices_raw; i++) {
        if (new_id[i] == 0) continue;
        for (int nb : adj[i]) {
            new_adj[new_id[i]].push_back(new_id[nb]);
        }
    }

    // count edges (each undirected edge counted once)
    long long num_edges = 0;
    for (int i = 1; i <= num_vertices; i++) {
        num_edges += (long long)new_adj[i].size();
    }
    num_edges /= 2;

    // Step 4: write id_map file (new_id -> original name)
    {
        // build reverse map: old_id -> name
        std::vector<std::string> id_to_name(num_vertices_raw + 1);
        for (const auto& kv : name_to_id) {
            id_to_name[kv.second] = kv.first;
        }
        std::ofstream fmap(map_file);
        if (!fmap) {
            std::cerr << "Error: cannot write map file: " << map_file << "\n";
            return 1;
        }
        fmap << "new_id\toriginal_name\n";
        for (int i = 1; i <= num_vertices_raw; i++) {
            if (new_id[i] == 0) continue;
            fmap << new_id[i] << "\t" << id_to_name[i] << "\n";
        }
    }

    // Step 5: write METIS format
    std::ofstream fout(output_file);
    if (!fout) {
        std::cerr << "Error: cannot write output file: " << output_file << "\n";
        return 1;
    }

    fout << num_vertices << " " << num_edges << "\n";
    for (int i = 1; i <= num_vertices; i++) {
        bool first = true;
        for (int nb : new_adj[i]) {
            if (!first) fout << " ";
            fout << nb;
            first = false;
        }
        fout << "\n";
    }

    int isolated = num_vertices_raw - num_vertices;
    std::cerr << "Done: " << num_vertices << " vertices, " << num_edges << " edges"
              << " (" << (full_clique ? "clique" : "star") << " topology)"
              << (isolated > 0 ? ", " + std::to_string(isolated) + " isolated nodes skipped" : "")
              << " -> " << output_file << "\n"
              << "ID map: " << map_file << "\n";
    return 0;
}
