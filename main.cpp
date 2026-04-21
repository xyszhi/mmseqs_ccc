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
#include <cstring>
#include <cassert>
#include <chrono>
#include <iomanip>
#include <sys/stat.h>

// ---------------------------------------------------------------------------
// Timing helpers
// ---------------------------------------------------------------------------
using Clock     = std::chrono::steady_clock;
using TimePoint = Clock::time_point;

static TimePoint g_start; // program start time

// Return seconds elapsed since program start
static double elapsed_s() {
    return std::chrono::duration<double>(Clock::now() - g_start).count();
}

// Return seconds between two time points
static double diff_s(TimePoint a, TimePoint b) {
    return std::chrono::duration<double>(b - a).count();
}

// Format seconds as "Xh Ym Zs" or "Xm Zs" or "Z.Xs"
static std::string fmt_time(double s) {
    char buf[64];
    if (s >= 3600.0) {
        int h = (int)(s / 3600);
        int m = (int)((s - h * 3600) / 60);
        int sec = (int)(s - h * 3600 - m * 60);
        snprintf(buf, sizeof(buf), "%dh %02dm %02ds", h, m, sec);
    } else if (s >= 60.0) {
        int m = (int)(s / 60);
        double sec = s - m * 60;
        snprintf(buf, sizeof(buf), "%dm %04.1fs", m, sec);
    } else {
        snprintf(buf, sizeof(buf), "%.2fs", s);
    }
    return buf;
}

// Format large numbers with commas: 1234567 -> "1,234,567"
static std::string fmt_num(long long n) {
    std::string s = std::to_string(n);
    int pos = (int)s.size() - 3;
    while (pos > 0) { s.insert(pos, ","); pos -= 3; }
    return s;
}

// Print a timestamped log line: "[+12.34s] message"
static void log(const std::string& msg) {
    std::cout << "[+" << std::fixed << std::setprecision(2)
              << elapsed_s() << "s] " << msg << std::endl;
}

// ---------------------------------------------------------------------------
// Lock-free Union-Find using atomic CAS
// ---------------------------------------------------------------------------
struct DSU {
    std::vector<std::atomic<int>> parent;
    std::vector<std::atomic<int>> rnk;

    explicit DSU(int n) : parent(n), rnk(n) {
        for (int i = 0; i < n; i++) {
            parent[i].store(i, std::memory_order_relaxed);
            rnk[i].store(0, std::memory_order_relaxed);
        }
    }

    int find(int x) {
        while (true) {
            int p = parent[x].load(std::memory_order_relaxed);
            if (p == x) return x;
            int pp = parent[p].load(std::memory_order_relaxed);
            parent[x].compare_exchange_weak(p, pp, std::memory_order_relaxed);
            x = pp;
        }
    }

    void unite(int a, int b) {
        while (true) {
            a = find(a);
            b = find(b);
            if (a == b) return;

            int ra = rnk[a].load(std::memory_order_relaxed);
            int rb = rnk[b].load(std::memory_order_relaxed);
            if (ra < rb) std::swap(a, b);

            int expected_b = b;
            if (!parent[b].compare_exchange_strong(expected_b, a,
                    std::memory_order_acq_rel, std::memory_order_relaxed))
                continue;

            if (ra == rb)
                rnk[a].compare_exchange_strong(ra, ra + 1,
                    std::memory_order_relaxed, std::memory_order_relaxed);
            return;
        }
    }
};

// ---------------------------------------------------------------------------
// Fast ASCII integer parser
// ---------------------------------------------------------------------------
static inline const char* parse_int(const char* p, const char* end, int& out) {
    while (p < end && (*p == ' ' || *p == '\t')) ++p;
    if (p >= end || *p < '0' || *p > '9') return nullptr;
    int v = 0;
    while (p < end && *p >= '0' && *p <= '9') {
        v = v * 10 + (*p - '0');
        ++p;
    }
    out = v;
    return p;
}

struct IndexEntry {
    int       query_id;
    long long local_offset;
    long long length;
};

int main(int argc, char* argv[]) {
    g_start = Clock::now();

    std::string db_prefix   = "test/filtered_db";
    std::string lookup_path = "test/protein_db.lookup";
    int num_threads = 8;
    if (argc > 1) db_prefix   = argv[1];
    if (argc > 2) lookup_path = argv[2];
    if (argc > 3) num_threads = std::stoi(argv[3]);

    std::cout << "========================================" << std::endl;
    std::cout << "  mmseqs_ccc  Connected-Component Clustering" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "  DB prefix   : " << db_prefix   << std::endl;
    std::cout << "  Lookup file : " << lookup_path  << std::endl;
    std::cout << "  Threads     : " << num_threads  << std::endl;
    std::cout << "========================================" << std::endl;

    // -----------------------------------------------------------------------
    // Phase 1: Load lookup
    // -----------------------------------------------------------------------
    log("Phase 1/4  Loading lookup file...");
    auto t0 = Clock::now();

    int max_id = 0;
    std::vector<std::pair<int,std::string>> id_name_pairs;
    {
        std::ifstream lf(lookup_path);
        if (!lf.is_open()) {
            std::cerr << "Warning: cannot open lookup file: " << lookup_path << std::endl;
        } else {
            std::string lline;
            lline.reserve(256);
            while (std::getline(lf, lline)) {
                if (lline.empty()) continue;
                const char* p   = lline.data();
                const char* end = p + lline.size();
                int id;
                p = parse_int(p, end, id);
                if (!p) continue;
                while (p < end && (*p == ' ' || *p == '\t')) ++p;
                const char* name_start = p;
                while (p < end && *p != '\r' && *p != '\n') ++p;
                if (id + 1 > max_id) max_id = id + 1;
                id_name_pairs.emplace_back(id, std::string(name_start, p));
            }
        }
    }

    std::vector<std::string> id2name;
    if (max_id > 0) {
        id2name.resize(max_id);
        for (auto& kv : id_name_pairs)
            id2name[kv.first] = std::move(kv.second);
        id_name_pairs.clear();
        id_name_pairs.shrink_to_fit();
    }

    {
        double dt = diff_s(t0, Clock::now());
        std::cout << "  -> Loaded " << fmt_num((long long)id2name.size())
                  << " sequence slots, max_id=" << fmt_num(max_id)
                  << "  (" << fmt_time(dt) << ")" << std::endl;
    }

    // -----------------------------------------------------------------------
    // Phase 2: Discover split files and read index
    // -----------------------------------------------------------------------
    log("Phase 2/4  Scanning split files and reading index...");
    t0 = Clock::now();

    int num_splits = 0;
    {
        struct stat st;
        while (stat((db_prefix + "." + std::to_string(num_splits)).c_str(), &st) == 0)
            num_splits++;
    }

    std::vector<long long> cum_sizes(num_splits + 1, 0);
    long long total_db_bytes = 0;
    for (int i = 0; i < num_splits; i++) {
        struct stat st;
        stat((db_prefix + "." + std::to_string(i)).c_str(), &st);
        cum_sizes[i + 1] = cum_sizes[i] + (long long)st.st_size;
    }
    total_db_bytes = cum_sizes[num_splits];

    std::vector<std::vector<IndexEntry>> file_entries(num_splits);
    long long total_index_entries = 0;
    {
        std::ifstream idx(db_prefix + ".index");
        if (!idx.is_open()) {
            std::cerr << "Cannot open index file: " << db_prefix << ".index" << std::endl;
            return 1;
        }
        std::string line;
        line.reserve(128);
        while (std::getline(idx, line)) {
            if (line.empty()) continue;
            const char* p   = line.data();
            const char* end = p + line.size();
            int query_id;
            p = parse_int(p, end, query_id);
            if (!p) continue;

            while (p < end && (*p == ' ' || *p == '\t')) ++p;
            long long offset = 0;
            while (p < end && *p >= '0' && *p <= '9') { offset = offset*10+(*p-'0'); ++p; }

            while (p < end && (*p == ' ' || *p == '\t')) ++p;
            long long length = 0;
            while (p < end && *p >= '0' && *p <= '9') { length = length*10+(*p-'0'); ++p; }

            int fid = (int)(std::upper_bound(cum_sizes.begin(), cum_sizes.end(), offset)
                            - cum_sizes.begin()) - 1;
            if (fid < 0) fid = 0;
            if (fid >= num_splits) fid = num_splits - 1;

            file_entries[fid].push_back({query_id, offset - cum_sizes[fid], length});
            total_index_entries++;
        }
    }

    for (int i = 0; i < num_splits; i++) {
        std::sort(file_entries[i].begin(), file_entries[i].end(),
                  [](const IndexEntry& a, const IndexEntry& b) {
                      return a.local_offset < b.local_offset;
                  });
    }

    {
        double dt = diff_s(t0, Clock::now());
        double db_gb = total_db_bytes / 1073741824.0;
        std::cout << "  -> " << num_splits << " split files, total size "
                  << std::fixed << std::setprecision(2) << db_gb << " GB"
                  << ";  " << fmt_num(total_index_entries) << " index entries"
                  << "  (" << fmt_time(dt) << ")" << std::endl;
    }

    // -----------------------------------------------------------------------
    // Phase 3: Parallel clustering (DSU)
    // -----------------------------------------------------------------------
    log("Phase 3/4  Building DSU with " + std::to_string(num_threads) + " threads...");
    t0 = Clock::now();

    if (max_id == 0) max_id = 20000000;
    DSU dsu(max_id);
    std::atomic<long long> total_edges{0};
    std::atomic<int>       splits_done{0};
    std::atomic<int>       next_file{0};

    // Progress reporter thread
    std::atomic<bool> progress_stop{false};
    std::thread progress_thread([&]() {
        while (!progress_stop.load()) {
            std::this_thread::sleep_for(std::chrono::seconds(30));
            if (progress_stop.load()) break;
            int done  = splits_done.load();
            long long edges = total_edges.load();
            double pct = num_splits > 0 ? 100.0 * done / num_splits : 0.0;
            double dt  = diff_s(t0, Clock::now());
            double eps = dt > 0 ? edges / dt : 0;
            std::cout << "  [progress] splits " << done << "/" << num_splits
                      << " (" << std::fixed << std::setprecision(1) << pct << "%)"
                      << "  edges=" << fmt_num(edges)
                      << "  speed=" << fmt_num((long long)eps) << " edges/s"
                      << "  elapsed=" << fmt_time(dt)
                      << std::endl;
        }
    });

    const int IO_BUF = 64 << 20; // 64 MB per thread

    auto worker = [&]() {
        std::vector<char> io_buf(IO_BUF);
        int fid;
        while ((fid = next_file.fetch_add(1)) < num_splits) {
            auto& entries = file_entries[fid];
            if (entries.empty()) {
                splits_done.fetch_add(1);
                continue;
            }

            std::ifstream sf(db_prefix + "." + std::to_string(fid), std::ios::binary);
            if (!sf.is_open()) {
                std::cerr << "Cannot open split file " << fid << std::endl;
                splits_done.fetch_add(1);
                continue;
            }

            long long local_edges = 0;
            long long buf_start = -1, buf_end = -1;

            auto ensure_buf = [&](long long need_start, long long need_len) -> bool {
                if (need_start >= buf_start && need_start + need_len <= buf_end)
                    return true;
                sf.seekg(need_start);
                if (!sf) return false;
                long long to_read = std::max(need_len, (long long)IO_BUF);
                sf.read(io_buf.data(), to_read);
                long long got = sf.gcount();
                if (got <= 0) return false;
                buf_start = need_start;
                buf_end   = need_start + got;
                return (got >= need_len);
            };

            for (auto& e : entries) {
                const char* block_ptr = nullptr;
                std::string fallback_buf;

                if (e.length <= IO_BUF) {
                    if (!ensure_buf(e.local_offset, e.length)) continue;
                    block_ptr = io_buf.data() + (e.local_offset - buf_start);
                } else {
                    fallback_buf.resize(e.length);
                    sf.seekg(e.local_offset);
                    sf.read(&fallback_buf[0], e.length);
                    if (sf.gcount() < e.length) continue;
                    block_ptr = fallback_buf.data();
                    buf_start = buf_end = -1;
                }

                const char* p   = block_ptr;
                const char* end = block_ptr + e.length;

                while (p < end) {
                    if (*p == '\0') { ++p; continue; }
                    const char* lend = p;
                    while (lend < end && *lend != '\n') ++lend;

                    int target_id;
                    const char* r = parse_int(p, lend, target_id);
                    if (r && target_id < max_id) {
                        dsu.unite(e.query_id, target_id);
                        local_edges++;
                    }
                    p = (lend < end) ? lend + 1 : end;
                }
            }

            total_edges += local_edges;
            int done = splits_done.fetch_add(1) + 1;
            double pct = 100.0 * done / num_splits;
            double dt  = diff_s(t0, Clock::now());
            double eps = dt > 0 ? total_edges.load() / dt : 0;

            // Print per-split line with progress info
            std::cout << "  split " << std::setw(4) << fid
                      << "  queries=" << std::setw(7) << entries.size()
                      << "  edges=" << std::setw(9) << fmt_num(local_edges)
                      << "  [" << std::setw(5) << std::fixed << std::setprecision(1) << pct << "%"
                      << "  " << done << "/" << num_splits << "]"
                      << "  speed=" << fmt_num((long long)eps) << " e/s"
                      << "  elapsed=" << fmt_time(dt)
                      << std::endl;

            // Free index entries immediately
            std::vector<IndexEntry>().swap(file_entries[fid]);
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(num_threads);
    for (int i = 0; i < num_threads; i++)
        threads.emplace_back(worker);
    for (auto& t : threads) t.join();

    progress_stop.store(true);
    progress_thread.join();

    {
        double dt = diff_s(t0, Clock::now());
        long long edges = total_edges.load();
        std::cout << "  -> Clustering done: " << fmt_num(edges) << " edges processed"
                  << "  avg speed=" << fmt_num((long long)(dt > 0 ? edges / dt : 0)) << " edges/s"
                  << "  (" << fmt_time(dt) << ")" << std::endl;
    }

    // -----------------------------------------------------------------------
    // Phase 4: Collect components and write output
    // -----------------------------------------------------------------------
    log("Phase 4/4  Collecting connected components...");
    t0 = Clock::now();

    std::vector<std::pair<int,int>> root_member;
    root_member.reserve(max_id);
    for (int i = 0; i < max_id; i++) {
        if (!id2name[i].empty()) {
            root_member.emplace_back(dsu.find(i), i);
        }
    }
    std::cout << "  -> Built " << fmt_num((long long)root_member.size())
              << " (root,member) pairs, sorting..." << std::endl;

    std::sort(root_member.begin(), root_member.end());

    int num_clusters = 0, singleton = 0, largest = 0;
    {
        int i = 0, n = (int)root_member.size();
        while (i < n) {
            int root = root_member[i].first;
            int j = i;
            while (j < n && root_member[j].first == root) ++j;
            int sz = j - i;
            num_clusters++;
            if (sz == 1) singleton++;
            if (sz > largest) largest = sz;
            i = j;
        }
    }

    std::cout << "  -> Clusters      : " << fmt_num(num_clusters) << std::endl;
    std::cout << "  -> Largest cluster: " << fmt_num(largest)      << std::endl;
    std::cout << "  -> Singletons    : " << fmt_num(singleton)     << std::endl;

    std::string out_path = db_prefix + "_clusters.tsv";
    std::cout << "  -> Writing output to: " << out_path << " ..." << std::endl;

    std::ofstream out(out_path);
    const int OUT_BUF = 8 << 20;
    std::vector<char> out_buf_storage(OUT_BUF);
    out.rdbuf()->pubsetbuf(out_buf_storage.data(), OUT_BUF);

    {
        int i = 0, n = (int)root_member.size();
        while (i < n) {
            int root = root_member[i].first;
            int j = i;
            while (j < n && root_member[j].first == root) ++j;
            const std::string& rep_name = id2name[root_member[i].second];
            for (int k = i; k < j; k++)
                out << rep_name << '\t' << id2name[root_member[k].second] << '\n';
            i = j;
        }
    }
    out.close();

    {
        double dt = diff_s(t0, Clock::now());
        std::cout << "  -> Output written  (" << fmt_time(dt) << ")" << std::endl;
    }

    std::cout << "========================================" << std::endl;
    log("All done.  Total wall time: " + fmt_time(elapsed_s()));
    std::cout << "========================================" << std::endl;

    return 0;
}
