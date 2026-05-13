#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <numeric>
#include <algorithm>
#include <sys/stat.h>
#include <cerrno>
#include <getopt.h>
#include <thread>
#include <mutex>
#include <atomic>

/**
 * Lock-free Union-Find using atomic CAS (from main.cpp)
 */
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

// 从行中解析 qn 和 tn，去除尾部 \r，返回 false 表示格式无效
static bool parse_edge(const std::string& line, const std::string& sep,
                        std::string& qn, std::string& tn) {
    if (line.empty()) return false;
    size_t pos = line.find(sep);
    if (pos == std::string::npos) return false;
    qn = line.substr(0, pos);
    size_t pos2 = line.find(sep, pos + sep.size());
    tn = (pos2 == std::string::npos)
         ? line.substr(pos + sep.size())
         : line.substr(pos + sep.size(), pos2 - (pos + sep.size()));
    // 去除 Windows 换行符 \r
    if (!qn.empty() && qn.back() == '\r') qn.pop_back();
    if (!tn.empty() && tn.back() == '\r') tn.pop_back();
    return !qn.empty() && !tn.empty();
}

void print_usage(const char* prog) {
    std::cerr << "Usage: " << prog << " --input <tsv_file> --output <out_dir> [--sep <sep>] [--threads <n>]\n"
              << "Options:\n"
              << "  -i, --input <file>    Input TSV file\n"
              << "  -o, --output <dir>    Output directory\n"
              << "  -s, --sep <char>      Separator (default: \\t)\n"
              << "  -t, --threads <int>   Number of threads (default: hardware_concurrency)\n"
              << "  -h, --help            Show this help message\n";
}

int main(int argc, char* argv[]) {
    std::string input_path;
    std::string output_dir;
    std::string sep = "\t";
    int num_threads = std::thread::hardware_concurrency();
    if (num_threads == 0) num_threads = 8;

    static struct option long_options[] = {
        {"input",   required_argument, 0, 'i'},
        {"output",  required_argument, 0, 'o'},
        {"sep",     required_argument, 0, 's'},
        {"threads", required_argument, 0, 't'},
        {"help",    no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "i:o:s:t:h", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'i': input_path = optarg; break;
            case 'o': output_dir = optarg; break;
            case 's': sep = optarg; if (sep == "\\t") sep = "\t"; break;
            case 't': try { num_threads = std::stoi(optarg); }
                      catch (...) { std::cerr << "Error: invalid value for --threads\n"; return 1; }
                      if (num_threads < 1) { std::cerr << "Error: --threads must be >= 1\n"; return 1; }
                      break;
            case 'h': print_usage(argv[0]); return 0;
            default: print_usage(argv[0]); return 1;
        }
    }

    if (input_path.empty() || output_dir.empty()) {
        print_usage(argv[0]);
        return 1;
    }

    if (sep.empty()) {
        std::cerr << "Error: separator cannot be empty." << std::endl;
        return 1;
    }

    if (mkdir(output_dir.c_str(), 0755) != 0 && errno != EEXIST) {
        std::cerr << "Error: cannot create output directory: " << output_dir << std::endl;
        return 1;
    }

    // -------------------------------------------------------------------------
    // Pass 1: 串行发现所有唯一节点，建立 name -> id 映射
    // -------------------------------------------------------------------------
    std::cerr << "Pass 1/2: Discovering unique nodes..." << std::endl;
    std::unordered_map<std::string, int> name2id;
    {
        std::ifstream infile(input_path);
        if (!infile.is_open()) {
            std::cerr << "Error: cannot open " << input_path << std::endl;
            return 1;
        }
        std::string line;
        int node_count = 0;
        while (std::getline(infile, line)) {
            std::string qn, tn;
            if (!parse_edge(line, sep, qn, tn)) continue;
            if (name2id.find(qn) == name2id.end()) name2id[qn] = node_count++;
            if (name2id.find(tn) == name2id.end()) name2id[tn] = node_count++;
        }
    }
    int total_nodes = (int)name2id.size();
    std::cerr << "Found " << total_nodes << " unique nodes." << std::endl;

    if (total_nodes == 0) {
        std::cerr << "Error: no valid edges found in input file." << std::endl;
        return 1;
    }

    // -------------------------------------------------------------------------
    // Pass 2 (DSU): 并行读取边，构建并查集
    // name2id 在此阶段只读，并发安全；使用 find() 而非 operator[] 避免静默插入
    // -------------------------------------------------------------------------
    std::cerr << "Building DSU (Parallel)..." << std::endl;
    DSU dsu(total_nodes);
    {
        std::ifstream infile(input_path);
        if (!infile.is_open()) {
            std::cerr << "Error: cannot open " << input_path << std::endl;
            return 1;
        }
        std::mutex file_mtx;
        auto worker = [&]() {
            std::string line;
            std::vector<std::string> local_lines;
            local_lines.reserve(1000);
            while (true) {
                {
                    std::lock_guard<std::mutex> lock(file_mtx);
                    for (int i = 0; i < 1000 && std::getline(infile, line); ++i)
                        local_lines.push_back(line);
                }
                if (local_lines.empty()) break;
                for (const auto& l : local_lines) {
                    std::string qn, tn;
                    if (!parse_edge(l, sep, qn, tn)) continue;
                    // 使用 find() 而非 operator[]，避免并发场景下的静默插入
                    auto it_q = name2id.find(qn);
                    auto it_t = name2id.find(tn);
                    if (it_q == name2id.end() || it_t == name2id.end()) continue;
                    dsu.unite(it_q->second, it_t->second);
                }
                local_lines.clear();
            }
        };
        std::vector<std::thread> threads;
        for (int i = 0; i < num_threads; ++i) threads.emplace_back(worker);
        for (auto& t : threads) t.join();
    }

    // -------------------------------------------------------------------------
    // 路径压缩：确保所有节点的根已完全收敛，再建立 root -> comp 映射
    // -------------------------------------------------------------------------
    // 先对所有节点做一次 find，使路径完全压缩
    for (int i = 0; i < total_nodes; ++i) dsu.find(i);

    std::unordered_map<int, int> root2comp;
    root2comp.reserve(total_nodes);
    int comp_count = 0;
    for (int i = 0; i < total_nodes; ++i) {
        int r = dsu.find(i);
        if (root2comp.find(r) == root2comp.end())
            root2comp[r] = comp_count++;
    }
    std::cerr << "Found " << comp_count << " connected components." << std::endl;

    // 预计算每个节点的 component id，避免 Pass 3 中重复查找并消除 operator[] 风险
    // node_comp[i] = component id of node i
    std::vector<int> node_comp(total_nodes);
    for (int i = 0; i < total_nodes; ++i) {
        node_comp[i] = root2comp[dsu.find(i)];
    }

    // -------------------------------------------------------------------------
    // Pass 3: 并行分发数据到各 component 文件
    // -------------------------------------------------------------------------
    std::cerr << "Pass 2/2: Distributing data..." << std::endl;
    {
        std::ifstream infile(input_path);
        if (!infile.is_open()) {
            std::cerr << "Error: cannot open " << input_path << std::endl;
            return 1;
        }
        std::mutex file_mtx;
        std::vector<std::mutex> comp_locks(comp_count);

        auto worker = [&]() {
            std::string line;
            std::vector<std::string> local_lines;
            local_lines.reserve(1000);
            std::unordered_map<int, std::string> thread_buffers;
            size_t buffered_bytes = 0;

            auto flush_buffers = [&]() {
                for (auto& kv : thread_buffers) {
                    int cid = kv.first;
                    std::lock_guard<std::mutex> lock(comp_locks[cid]);
                    std::string out_file = output_dir + "/component_" + std::to_string(cid) + ".tsv";
                    std::ofstream out(out_file, std::ios::app);
                    if (!out.is_open()) {
                        std::cerr << "Error: cannot open output file: " << out_file << std::endl;
                    } else {
                        out << kv.second;
                        if (out.fail())
                            std::cerr << "Error: write failed for: " << out_file << std::endl;
                    }
                }
                thread_buffers.clear();
                buffered_bytes = 0;
            };

            while (true) {
                {
                    std::lock_guard<std::mutex> lock(file_mtx);
                    for (int i = 0; i < 1000 && std::getline(infile, line); ++i)
                        local_lines.push_back(line);
                }
                if (local_lines.empty()) break;
                for (const auto& l : local_lines) {
                    std::string qn, tn;
                    if (!parse_edge(l, sep, qn, tn)) continue;
                    // 使用 find() 安全查找，跳过不在 name2id 中的节点
                    auto it = name2id.find(qn);
                    if (it == name2id.end()) continue;
                    int cid = node_comp[it->second];
                    // 写出原始行（去掉可能的尾部 \r 后重建，保证输出干净）
                    thread_buffers[cid] += l;
                    // 去除尾部 \r，统一用 \n 结尾
                    if (!thread_buffers[cid].empty() && thread_buffers[cid].back() == '\r')
                        thread_buffers[cid].back() = '\n';
                    else
                        thread_buffers[cid] += '\n';
                    buffered_bytes += l.size() + 1;
                    if (buffered_bytes > 10 * 1024 * 1024) flush_buffers();
                }
                local_lines.clear();
            }
            flush_buffers();
        };

        std::vector<std::thread> threads;
        for (int i = 0; i < num_threads; ++i) threads.emplace_back(worker);
        for (auto& t : threads) t.join();
    }

    std::cerr << "Done." << std::endl;
    return 0;
}
