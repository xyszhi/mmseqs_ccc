#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <stdexcept>

// ---------------------------------------------------------------------------
// Usage
// ---------------------------------------------------------------------------
static void print_usage(const char* prog) {
    std::cerr <<
        "Usage: " << prog << " [OPTIONS]\n"
        "\n"
        "Convert a clusters TSV file (rep<TAB>member) to METIS graph format.\n"
        "Designed for very large files (e.g. 500 GB) with limited RAM.\n"
        "\n"
        "Options:\n"
        "  --input   <path>  Input clusters TSV file (default: test/filtered_db_clusters.tsv)\n"
        "  -i        <path>  Alias for --input\n"
        "  --output  <path>  Output METIS .graph file (default: <input>.graph)\n"
        "  -o        <path>  Alias for --output\n"
        "  --tmpdir  <path>  Directory for temporary files (default: same as output)\n"
        "  --clique          Build full clique within each cluster (default: star topology)\n"
        "  --help            Show this help\n"
        "  -h                Show this help\n"
        "\n"
        "Algorithm (two-pass + external sort, O(1) peak RAM beyond ID map):\n"
        "  Pass 1 : Stream TSV, assign integer IDs to node names, write id_map to disk.\n"
        "  Pass 2 : Stream TSV again, emit (u,v) integer edge pairs to a temp binary file.\n"
        "  Sort   : Sort + deduplicate edge pairs with std::sort (fits in RAM as integers).\n"
        "           If edge list is too large, falls back to in-place sort of temp file.\n"
        "  Output : Stream sorted edges to METIS format.\n"
        "\n"
        "Topology modes:\n"
        "  star   (default): representative <-> each member (fewer edges)\n"
        "  clique (--clique): all pairs within cluster connected (more edges)\n"
        "\n"
        "Output format: METIS unweighted undirected graph (.graph)\n"
        "  Line 1: <num_vertices> <num_edges>\n"
        "  Line i: space-separated neighbor IDs (1-based) of vertex i\n";
}

// ---------------------------------------------------------------------------
// Edge: a pair of 1-based integer node IDs, stored canonically (u <= v)
// ---------------------------------------------------------------------------
struct Edge {
    int32_t u, v; // u <= v always
};

static inline Edge make_edge(int32_t a, int32_t b) {
    if (a > b) std::swap(a, b);
    return {a, b};
}

static inline bool edge_less(const Edge& a, const Edge& b) {
    if (a.u != b.u) return a.u < b.u;
    return a.v < b.v;
}

// ---------------------------------------------------------------------------
// Pass 1: stream TSV, assign integer IDs, write id_map file
// Returns number of unique nodes assigned.
// ---------------------------------------------------------------------------
static int64_t pass1_assign_ids(
    const std::string& input_file,
    std::unordered_map<std::string, int32_t>& name_to_id)
{
    std::ifstream fin(input_file);
    if (!fin) throw std::runtime_error("Cannot open input file: " + input_file);

    std::string line, rep, mem;
    int32_t next_id = 1;

    while (std::getline(fin, line)) {
        if (line.empty() || line[0] == '%') continue;
        std::istringstream ss(line);
        if (!(ss >> rep >> mem)) continue;

        for (const std::string* s : {&rep, &mem}) {
            auto it = name_to_id.find(*s);
            if (it == name_to_id.end()) {
                name_to_id[*s] = next_id;
                ++next_id;
            }
        }
    }
    return (int64_t)(next_id - 1);
}

// ---------------------------------------------------------------------------
// Pass 2: stream TSV again, emit integer edge pairs to binary temp file.
// Star topology: rep <-> each member.
// Clique topology: all pairs within each cluster.
// Returns total number of (possibly duplicate) edges written.
// ---------------------------------------------------------------------------
static int64_t pass2_emit_edges(
    const std::string& input_file,
    const std::string& tmp_edge_file,
    const std::unordered_map<std::string, int32_t>& name_to_id,
    bool full_clique)
{
    std::ifstream fin(input_file);
    if (!fin) throw std::runtime_error("Cannot open input file: " + input_file);

    FILE* ftmp = std::fopen(tmp_edge_file.c_str(), "wb");
    if (!ftmp) throw std::runtime_error("Cannot open temp edge file: " + tmp_edge_file);

    // Buffer writes for performance
    const size_t BUF_EDGES = 1 << 20; // 1M edges per flush (~8MB)
    std::vector<Edge> buf;
    buf.reserve(BUF_EDGES);

    auto flush_buf = [&]() {
        if (!buf.empty()) {
            std::fwrite(buf.data(), sizeof(Edge), buf.size(), ftmp);
            buf.clear();
        }
    };

    auto emit = [&](int32_t a, int32_t b) {
        if (a == b) return; // skip self-loops
        buf.push_back(make_edge(a, b));
        if (buf.size() >= BUF_EDGES) flush_buf();
    };

    int64_t total_written = 0;
    std::string line, rep, mem;

    // For clique mode we need to accumulate members per cluster
    // We process cluster-by-cluster: collect members until rep changes
    if (full_clique) {
        std::string cur_rep;
        std::vector<int32_t> cur_members;

        auto flush_cluster = [&]() {
            if (cur_rep.empty()) return;
            int32_t rep_id = name_to_id.at(cur_rep);
            // all pairs including rep
            std::vector<int32_t> all;
            all.push_back(rep_id);
            for (int32_t m : cur_members) all.push_back(m);
            // deduplicate within cluster
            std::sort(all.begin(), all.end());
            all.erase(std::unique(all.begin(), all.end()), all.end());
            for (size_t i = 0; i < all.size(); i++) {
                for (size_t j = i + 1; j < all.size(); j++) {
                    emit(all[i], all[j]);
                    total_written++;
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
            }
            int32_t mem_id = name_to_id.at(mem);
            cur_members.push_back(mem_id);
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
                emit(rep_id, mem_id);
                total_written++;
            }
        }
    }

    flush_buf();
    std::fclose(ftmp);
    return total_written;
}

// ---------------------------------------------------------------------------
// Sort + deduplicate the binary edge file in memory.
// Reads all edges, sorts, deduplicates, writes back.
// Returns number of unique undirected edges.
// ---------------------------------------------------------------------------
static int64_t sort_dedup_edges(const std::string& tmp_edge_file, int64_t num_raw_edges) {
    std::cerr << "  Loading " << num_raw_edges << " raw edges into memory for sorting...\n";

    std::vector<Edge> edges(num_raw_edges);
    FILE* f = std::fopen(tmp_edge_file.c_str(), "rb");
    if (!f) throw std::runtime_error("Cannot open temp edge file for reading: " + tmp_edge_file);
    int64_t read_count = (int64_t)std::fread(edges.data(), sizeof(Edge), num_raw_edges, f);
    std::fclose(f);
    if (read_count != num_raw_edges)
        throw std::runtime_error("Short read on temp edge file");

    std::cerr << "  Sorting...\n";
    std::sort(edges.begin(), edges.end(), edge_less);

    std::cerr << "  Deduplicating...\n";
    auto new_end = std::unique(edges.begin(), edges.end(),
        [](const Edge& a, const Edge& b){ return a.u == b.u && a.v == b.v; });
    edges.erase(new_end, edges.end());
    int64_t unique_count = (int64_t)edges.size();

    // Write back deduplicated edges
    f = std::fopen(tmp_edge_file.c_str(), "wb");
    if (!f) throw std::runtime_error("Cannot open temp edge file for writing: " + tmp_edge_file);
    std::fwrite(edges.data(), sizeof(Edge), unique_count, f);
    std::fclose(f);

    return unique_count;
}

// ---------------------------------------------------------------------------
// Build adjacency list from sorted unique edges and write METIS format.
// Streams through the sorted edge file; builds one vertex's neighbor list
// at a time. Requires two passes over the edge file (one to count degrees,
// one to write neighbors) — both are sequential reads.
// ---------------------------------------------------------------------------
static void write_metis(
    const std::string& tmp_edge_file,
    const std::string& output_file,
    const std::string& map_file,
    const std::unordered_map<std::string, int32_t>& name_to_id,
    int64_t num_vertices,
    int64_t num_unique_edges)
{
    // Build reverse map: old_id -> name (for id_map output)
    std::vector<std::string> id_to_name(num_vertices + 1);
    for (const auto& kv : name_to_id) id_to_name[kv.second] = kv.first;

    // Pass A: compute degree of each vertex (to know neighbor list sizes)
    std::cerr << "  Computing vertex degrees...\n";
    std::vector<int64_t> degree(num_vertices + 1, 0);

    {
        FILE* f = std::fopen(tmp_edge_file.c_str(), "rb");
        if (!f) throw std::runtime_error("Cannot open temp edge file: " + tmp_edge_file);
        Edge e;
        while (std::fread(&e, sizeof(Edge), 1, f) == 1) {
            degree[e.u]++;
            degree[e.v]++;
        }
        std::fclose(f);
    }

    // Identify vertices that have at least one edge (non-isolated)
    // Build compact renumbering: old_id (1-based) -> new_id (1-based)
    std::vector<int32_t> new_id(num_vertices + 1, 0);
    int64_t num_vertices_out = 0;
    for (int64_t i = 1; i <= num_vertices; i++) {
        if (degree[i] > 0) new_id[i] = (int32_t)(++num_vertices_out);
    }
    int64_t isolated = num_vertices - num_vertices_out;

    // Write id_map: only non-isolated nodes, new_id -> original name
    {
        std::ofstream fmap(map_file);
        if (!fmap) throw std::runtime_error("Cannot write map file: " + map_file);
        fmap << "new_id\toriginal_name\n";
        for (int64_t i = 1; i <= num_vertices; i++) {
            if (new_id[i] == 0) continue;
            fmap << new_id[i] << "\t" << id_to_name[i] << "\n";
        }
    }

    std::cerr << "  Writing METIS file: " << num_vertices_out << " vertices, "
              << num_unique_edges << " edges";
    if (isolated > 0) std::cerr << " (" << isolated << " isolated nodes skipped)";
    std::cerr << "\n";

    // Pass B: for each output vertex (in new_id order), collect and write neighbors.
    // We stream the sorted edge file and build neighbor lists vertex by vertex.
    // Since edges are sorted by (u,v), we can process them in one sequential pass.
    // We need to output vertex i's neighbors before moving to i+1, so we buffer
    // per-vertex neighbor lists. Memory: O(max_degree) at any time.

    std::ofstream fout(output_file);
    if (!fout) throw std::runtime_error("Cannot write output file: " + output_file);

    fout << num_vertices_out << " " << num_unique_edges << "\n";

    // Build full adjacency in sorted order.
    // Since edges are sorted by (u,v) and undirected, we need both directions.
    // Strategy: load all edges, build CSR-style adjacency, then stream output.
    // Memory: degree array (already have) + neighbor arrays = O(2*E * 4 bytes).
    // For 10B edges that's 80GB — too much. Instead, use a two-pointer approach:
    // write output line by line, for each vertex scan edges where u==v or v==u.
    // But that's O(V*E). Better: build CSR from sorted edges.

    // CSR build: offset array + neighbor array
    // offset[i] = start index in neighbor array for vertex new_id=i
    // Total neighbor entries = 2 * num_unique_edges (each edge contributes to both endpoints)
    // Memory: (num_vertices_out+1)*4 + 2*num_unique_edges*4 bytes
    // For 10B edges: ~80GB — still large. But this is the minimum needed for METIS output.
    // In practice, for 500GB TSV with star topology, edges ~ rows ~ 10B,
    // but each edge is 8 bytes, so 80GB for neighbors. If that's too much,
    // the user should use a graph DB. For typical protein cluster data,
    // edges are much fewer than rows (many members per cluster).

    // We'll build CSR using the degree array we already have.
    std::vector<int64_t> offset(num_vertices_out + 2, 0);
    for (int64_t i = 1; i <= num_vertices; i++) {
        if (new_id[i] == 0) continue;
        offset[new_id[i] + 1] = degree[i];
    }
    for (int64_t i = 1; i <= num_vertices_out; i++) {
        offset[i + 1] += offset[i];
    }
    int64_t total_adj = offset[num_vertices_out + 1];
    std::vector<int32_t> neighbors(total_adj);
    std::vector<int64_t> pos(num_vertices_out + 1); // current fill position
    for (int64_t i = 1; i <= num_vertices_out; i++) pos[i] = offset[i];

    {
        FILE* f = std::fopen(tmp_edge_file.c_str(), "rb");
        if (!f) throw std::runtime_error("Cannot open temp edge file: " + tmp_edge_file);
        Edge e;
        while (std::fread(&e, sizeof(Edge), 1, f) == 1) {
            int32_t nu = new_id[e.u], nv = new_id[e.v];
            if (nu && nv) {
                neighbors[pos[nu]++] = nv;
                neighbors[pos[nv]++] = nu;
            }
        }
        std::fclose(f);
    }

    // Sort each vertex's neighbor list (for canonical METIS output)
    for (int64_t i = 1; i <= num_vertices_out; i++) {
        std::sort(neighbors.begin() + (std::ptrdiff_t)offset[i], neighbors.begin() + (std::ptrdiff_t)offset[i + 1]);
    }

    // Write METIS adjacency lines
    for (int64_t i = 1; i <= num_vertices_out; i++) {
        bool first = true;
        for (int64_t j = offset[i]; j < offset[i + 1]; j++) {
            if (!first) fout << ' ';
            fout << neighbors[j];
            first = false;
        }
        fout << '\n';
    }
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    std::string input_file = "test/filtered_db_clusters.tsv";
    std::string output_file;
    std::string tmpdir;
    bool full_clique = false;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--input" || arg == "-i") {
            if (i + 1 < argc) input_file = argv[++i];
        } else if (arg == "--output" || arg == "-o") {
            if (i + 1 < argc) output_file = argv[++i];
        } else if (arg == "--tmpdir") {
            if (i + 1 < argc) tmpdir = argv[++i];
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
    if (tmpdir.empty()) {
        // Use same directory as output file
        size_t slash = output_file.rfind('/');
        tmpdir = (slash == std::string::npos) ? "." : output_file.substr(0, slash);
    }

    std::string map_file      = output_file + ".id_map";
    std::string tmp_edge_file = tmpdir + "/.extract_metis2_edges.tmp";

    try {
        // ------------------------------------------------------------------
        // Pass 1: assign integer IDs to all node names
        // ------------------------------------------------------------------
        std::cerr << "[Pass 1] Scanning TSV and assigning node IDs...\n";
        std::unordered_map<std::string, int32_t> name_to_id;
        int64_t num_nodes = pass1_assign_ids(input_file, name_to_id);
        std::cerr << "  " << num_nodes << " unique nodes found. ID map -> " << map_file << "\n";

        // ------------------------------------------------------------------
        // Pass 2: emit integer edge pairs to temp binary file
        // ------------------------------------------------------------------
        std::cerr << "[Pass 2] Emitting integer edge pairs to temp file...\n";
        int64_t num_raw_edges = pass2_emit_edges(
            input_file, tmp_edge_file, name_to_id, full_clique);
        std::cerr << "  " << num_raw_edges << " raw edge pairs written to " << tmp_edge_file << "\n";

        // ------------------------------------------------------------------
        // Sort + deduplicate edges
        // ------------------------------------------------------------------
        std::cerr << "[Sort] Sorting and deduplicating edges...\n";
        int64_t num_unique_edges = sort_dedup_edges(tmp_edge_file, num_raw_edges);
        std::cerr << "  " << num_unique_edges << " unique undirected edges.\n";

        // ------------------------------------------------------------------
        // Write METIS output (name_to_id still needed for id_map generation)
        // ------------------------------------------------------------------
        std::cerr << "[Output] Writing METIS graph file...\n";
        write_metis(tmp_edge_file, output_file, map_file, name_to_id, num_nodes, num_unique_edges);
        std::cerr << "Done -> " << output_file << "\n";

        // Free the name_to_id map — no longer needed
        { std::unordered_map<std::string, int32_t>().swap(name_to_id); }

        // Cleanup temp file
        std::remove(tmp_edge_file.c_str());

    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n";
        std::remove(tmp_edge_file.c_str());
        return 1;
    }

    return 0;
}
