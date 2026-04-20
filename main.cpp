#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <mutex>
#include <thread>
#include <atomic>
#include <cstdint>
#include <sys/stat.h>

// Thread-safe Union-Find with path compression and union by rank
struct DSU {
    std::vector<int> parent, rank_;

    explicit DSU(int n) : parent(n), rank_(n, 0) {
        for (int i = 0; i < n; i++) parent[i] = i;
    }

    int find(int x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]]; // path halving
            x = parent[x];
        }
        return x;
    }

    void unite(int a, int b) {
        a = find(a); b = find(b);
        if (a == b) return;
        if (rank_[a] < rank_[b]) std::swap(a, b);
        parent[b] = a;
        if (rank_[a] == rank_[b]) rank_[a]++;
    }
};

struct IndexEntry {
    int query_id;
    long long local_offset;
    long long length;
};

int main(int argc, char* argv[]) {
    std::string db_prefix = "test/filtered_db";
    std::string lookup_path = "test/protein_db.lookup";
    int num_threads = 8;
    if (argc > 1) db_prefix = argv[1];
    if (argc > 2) lookup_path = argv[2];
    if (argc > 3) num_threads = std::stoi(argv[3]);

    // Load ID -> sequence name mapping from lookup file
    std::unordered_map<int, std::string> id2name;
    {
        std::ifstream lf(lookup_path);
        if (lf.is_open()) {
            std::string lline;
            while (std::getline(lf, lline)) {
                if (lline.empty()) continue;
                std::istringstream ls(lline);
                int id; std::string name;
                ls >> id >> name;
                id2name[id] = name;
            }
            std::cout << "Loaded " << id2name.size() << " sequence names from lookup." << std::endl;
        } else {
            std::cerr << "Warning: cannot open lookup file: " << lookup_path << std::endl;
        }
    }

    int max_id = (int)id2name.size();

    // Count split files and compute cumulative sizes
    int num_splits = 0;
    {
        struct stat st;
        while (stat((db_prefix + "." + std::to_string(num_splits)).c_str(), &st) == 0) num_splits++;
    }
    std::cout << "Found " << num_splits << " split files." << std::endl;

    std::vector<long long> cum_sizes(num_splits + 1, 0);
    for (int i = 0; i < num_splits; i++) {
        struct stat st;
        stat((db_prefix + "." + std::to_string(i)).c_str(), &st);
        cum_sizes[i + 1] = cum_sizes[i] + (long long)st.st_size;
    }

    // Read index and group entries by split file, sorted by local offset
    std::vector<std::vector<IndexEntry>> file_entries(num_splits);
    {
        std::ifstream idx(db_prefix + ".index");
        if (!idx.is_open()) {
            std::cerr << "Cannot open index file: " << db_prefix << ".index" << std::endl;
            return 1;
        }
        std::string line;
        while (std::getline(idx, line)) {
            if (line.empty()) continue;
            std::istringstream ss(line);
            int query_id; long long offset, length;
            ss >> query_id >> offset >> length;

            // Binary search for file_id
            int fid = (int)(std::upper_bound(cum_sizes.begin(), cum_sizes.end(), offset) - cum_sizes.begin()) - 1;
            if (fid < 0) fid = 0;
            if (fid >= num_splits) fid = num_splits - 1;

            file_entries[fid].push_back({query_id, offset - cum_sizes[fid], length});
        }
    }

    // Sort each file's entries by local offset for sequential reads
    for (int i = 0; i < num_splits; i++) {
        std::sort(file_entries[i].begin(), file_entries[i].end(),
                  [](const IndexEntry& a, const IndexEntry& b) { return a.local_offset < b.local_offset; });
    }

    std::cout << "Index loaded and grouped. Starting clustering..." << std::endl;

    // DSU over all proteins
    if (max_id == 0) max_id = 20000000; // fallback if no lookup
    DSU dsu(max_id);
    std::mutex dsu_mutex;
    std::atomic<long long> total_edges{0};

    // Process split files in parallel
    std::vector<int> file_ids(num_splits);
    for (int i = 0; i < num_splits; i++) file_ids[i] = i;

    std::atomic<int> next_file{0};

    auto worker = [&]() {
        int fid;
        while ((fid = next_file.fetch_add(1)) < num_splits) {
            auto& entries = file_entries[fid];
            if (entries.empty()) continue;

            std::ifstream sf(db_prefix + "." + std::to_string(fid), std::ios::binary);
            if (!sf.is_open()) {
                std::cerr << "Cannot open split file " << fid << std::endl;
                continue;
            }

            long long local_edges = 0;
            // Collect edges locally to reduce lock contention
            std::vector<std::pair<int,int>> edges;
            edges.reserve(entries.size() * 4);

            for (auto& e : entries) {
                sf.seekg(e.local_offset);
                std::string block(e.length, '\0');
                sf.read(&block[0], e.length);

                std::istringstream block_ss(block);
                std::string hit_line;
                while (std::getline(block_ss, hit_line)) {
                    if (hit_line.empty() || hit_line[0] == '\0') continue;
                    while (!hit_line.empty() && (hit_line.back() == '\r' || hit_line.back() == '\0'))
                        hit_line.pop_back();
                    if (hit_line.empty()) continue;
                    std::istringstream hs(hit_line);
                    int target_id;
                    hs >> target_id;
                    if (target_id < max_id)
                        edges.push_back({e.query_id, target_id});
                    local_edges++;
                }
            }

            // Apply edges to DSU under lock
            {
                std::lock_guard<std::mutex> lock(dsu_mutex);
                for (auto& e2 : edges)
                    dsu.unite(e2.first, e2.second);
            }
            total_edges += local_edges;
            std::cout << "  Split " << fid << " done: " << entries.size() << " queries, " << local_edges << " edges." << std::endl;
        }
    };

    std::vector<std::thread> threads;
    for (int i = 0; i < num_threads; i++)
        threads.emplace_back(worker);
    for (auto& t : threads) t.join();

    std::cout << "Total edges processed: " << total_edges.load() << std::endl;

    // Collect connected components
    std::unordered_map<int, std::vector<int>> components;
    for (int i = 0; i < max_id; i++) {
        if (dsu.find(i) == i || dsu.parent[i] != i) {
            // only include nodes that appear in lookup
            if (id2name.count(i))
                components[dsu.find(i)].push_back(i);
        }
    }
    // Re-collect properly: iterate all known nodes
    components.clear();
    for (auto it = id2name.begin(); it != id2name.end(); ++it) {
        components[dsu.find(it->first)].push_back(it->first);
    }

    std::cout << "Number of connected components (clusters): " << components.size() << std::endl;

    // Statistics
    std::vector<int> sizes;
    sizes.reserve(components.size());
    for (auto it = components.begin(); it != components.end(); ++it)
        sizes.push_back((int)it->second.size());
    std::sort(sizes.rbegin(), sizes.rend());

    int singleton = 0;
    for (int s : sizes) if (s == 1) singleton++;
    std::cout << "Largest cluster size: " << sizes[0] << std::endl;
    std::cout << "Singleton clusters: " << singleton << std::endl;

    // Output clusters to file (with sequence names if available)
    std::string out_path = db_prefix + "_clusters.tsv";
    std::ofstream out(out_path);
    auto name_of = [&](int id) -> std::string {
        auto it = id2name.find(id);
        return (it != id2name.end()) ? it->second : std::to_string(id);
    };
    for (auto it = components.begin(); it != components.end(); ++it) {
        std::vector<int> sorted_members = it->second;
        std::sort(sorted_members.begin(), sorted_members.end());
        std::string root_name = name_of(it->first);
        for (int m : sorted_members)
            out << root_name << "\t" << name_of(m) << "\n";
    }
    out.close();
    std::cout << "Cluster results written to: " << out_path << std::endl;

    return 0;
}
