#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <map>
#include <cstdlib>

// ---------------------------------------------------------------------------
// Usage
// ---------------------------------------------------------------------------
static void print_usage(const char* prog) {
    std::cerr <<
        "Usage: " << prog << " [OPTIONS]\n"
        "\n"
        "Split a cluster alignment TSV file into per-partition sub-files\n"
        "based on KaMinPar partition results.\n"
        "\n"
        "Options:\n"
        "  --tsv       <path>  Input alignment TSV file (required)\n"
        "  --idmap     <path>  METIS id_map file (<graph>.metis.id_map) (required)\n"
        "  --partition <path>  Partition assignment file (partition.txt) (required)\n"
        "  --out-dir   <path>  Output directory (default: same as --tsv)\n"
        "  --help              Show this help\n"
        "\n"
        "Output files are named <tsv_basename>_K<block>.tsv in the output directory.\n"
        "Rows where query and target belong to different partitions are skipped.\n";
}

// ---------------------------------------------------------------------------
// Derive output directory from a file path (directory part)
// ---------------------------------------------------------------------------
static std::string dir_of(const std::string& path) {
    size_t slash = path.rfind('/');
    if (slash == std::string::npos) return ".";
    return path.substr(0, slash);
}

// ---------------------------------------------------------------------------
// Derive basename without extension from a file path
// ---------------------------------------------------------------------------
static std::string basename_noext(const std::string& path) {
    size_t slash = path.rfind('/');
    std::string base = (slash == std::string::npos) ? path : path.substr(slash + 1);
    size_t dot = base.rfind('.');
    if (dot != std::string::npos) base = base.substr(0, dot);
    return base;
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    std::string tsv_file;
    std::string idmap_file;
    std::string partition_file;
    std::string out_dir;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--tsv") {
            if (i + 1 < argc) tsv_file = argv[++i];
        } else if (arg == "--idmap") {
            if (i + 1 < argc) idmap_file = argv[++i];
        } else if (arg == "--partition") {
            if (i + 1 < argc) partition_file = argv[++i];
        } else if (arg == "--out-dir") {
            if (i + 1 < argc) out_dir = argv[++i];
        } else if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            return 0;
        } else {
            std::cerr << "Unknown argument: " << arg << "\n";
            print_usage(argv[0]);
            return 1;
        }
    }

    if (tsv_file.empty() || idmap_file.empty() || partition_file.empty()) {
        std::cerr << "Error: --tsv, --idmap and --partition are required.\n\n";
        print_usage(argv[0]);
        return 1;
    }

    if (out_dir.empty()) {
        out_dir = dir_of(tsv_file);
    }

    // -----------------------------------------------------------------------
    // Step 1: Load id_map: new_id (1-based) -> protein name
    // File format: header line, then "<new_id>\t<original_name>"
    // -----------------------------------------------------------------------
    std::cerr << "[Step 1] Loading id_map: " << idmap_file << "\n";
    // We need: protein_name -> new_id (for lookup)
    std::unordered_map<std::string, int32_t> name_to_nodeid;
    {
        std::ifstream fin(idmap_file);
        if (!fin) {
            std::cerr << "Error: cannot open id_map file: " << idmap_file << "\n";
            return 1;
        }
        std::string line;
        // skip header
        std::getline(fin, line);
        while (std::getline(fin, line)) {
            if (line.empty()) continue;
            std::istringstream ss(line);
            int32_t nid;
            std::string name;
            if (!(ss >> nid >> name)) continue;
            name_to_nodeid[name] = nid;
        }
    }
    std::cerr << "  " << name_to_nodeid.size() << " nodes loaded.\n";

    // -----------------------------------------------------------------------
    // Step 2: Load partition.txt: node_id (1-based) -> block_id
    // File format: comment lines starting with '#', then one block ID per line.
    // Line order (skipping comments) gives node_id starting from 1.
    // -----------------------------------------------------------------------
    std::cerr << "[Step 2] Loading partition: " << partition_file << "\n";
    // node_partition[i] = block for node (i+1), i.e. 0-based index
    std::vector<int32_t> node_partition;
    {
        std::ifstream fin(partition_file);
        if (!fin) {
            std::cerr << "Error: cannot open partition file: " << partition_file << "\n";
            return 1;
        }
        std::string line;
        while (std::getline(fin, line)) {
            if (line.empty() || line[0] == '#') continue;
            int32_t block = std::stoi(line);
            node_partition.push_back(block);
        }
    }
    std::cerr << "  " << node_partition.size() << " node assignments loaded.\n";

    // Helper: get block for a protein name; returns -1 if not found
    auto get_block = [&](const std::string& name) -> int32_t {
        auto it = name_to_nodeid.find(name);
        if (it == name_to_nodeid.end()) return -1;
        int32_t nid = it->second; // 1-based
        if (nid < 1 || nid > (int32_t)node_partition.size()) return -1;
        return node_partition[nid - 1];
    };

    // -----------------------------------------------------------------------
    // Step 3: Stream TSV, split rows into per-partition output files
    // -----------------------------------------------------------------------
    std::cerr << "[Step 3] Splitting TSV: " << tsv_file << "\n";

    std::ifstream fin(tsv_file);
    if (!fin) {
        std::cerr << "Error: cannot open TSV file: " << tsv_file << "\n";
        return 1;
    }

    std::string base = basename_noext(tsv_file);
    // Output files opened on demand
    std::map<int32_t, std::ofstream*> out_files;
    std::string header_line;

    auto get_outfile = [&](int32_t block) -> std::ofstream& {
        auto it = out_files.find(block);
        if (it != out_files.end()) return *(it->second);
        std::string path = out_dir + "/" + base + "_K" + std::to_string(block) + ".tsv";
        auto* ofs = new std::ofstream(path);
        if (!*ofs) {
            std::cerr << "Error: cannot open output file: " << path << "\n";
            std::exit(1);
        }
        // Write header if present
        if (!header_line.empty()) {
            *ofs << header_line << "\n";
        }
        out_files[block] = ofs;
        std::cerr << "  Opened output: " << path << "\n";
        return *ofs;
    };

    std::string line;
    long long total = 0, written = 0, skipped_cross = 0, skipped_unknown = 0;
    bool first_line = true;

    while (std::getline(fin, line)) {
        if (line.empty()) continue;

        // First line: treat as header if it starts with non-data character
        // or if first field doesn't look like a protein name with '|'
        // Simple heuristic: check if first token is "query"
        if (first_line) {
            first_line = false;
            std::istringstream ss(line);
            std::string tok;
            ss >> tok;
            if (tok == "query") {
                header_line = line;
                continue;
            }
            // Not a header, fall through to process as data
        }

        total++;

        // Extract query (col 1) and target (col 2)
        const char* p = line.c_str();
        const char* end = p + line.size();

        // find tab after col1
        const char* tab1 = (const char*)memchr(p, '\t', end - p);
        if (!tab1) { skipped_unknown++; continue; }
        std::string query(p, tab1);

        const char* col2_start = tab1 + 1;
        const char* tab2 = (const char*)memchr(col2_start, '\t', end - col2_start);
        std::string target = tab2 ? std::string(col2_start, tab2) : std::string(col2_start, end);

        int32_t bq = get_block(query);
        int32_t bt = get_block(target);

        if (bq < 0 || bt < 0) {
            skipped_unknown++;
            continue;
        }
        if (bq != bt) {
            skipped_cross++;
            continue;
        }

        get_outfile(bq) << line << "\n";
        written++;
    }

    // Close all output files
    for (auto& kv : out_files) {
        kv.second->close();
        delete kv.second;
    }

    std::cerr << "Done.\n";
    std::cerr << "  Total data rows : " << total << "\n";
    std::cerr << "  Written         : " << written << "\n";
    std::cerr << "  Skipped (cross) : " << skipped_cross << "\n";
    std::cerr << "  Skipped (unknown): " << skipped_unknown << "\n";
    std::cerr << "  Partitions output: " << out_files.size() << "\n";

    return 0;
}
