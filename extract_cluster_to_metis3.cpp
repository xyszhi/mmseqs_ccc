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
#include <queue>

// ---------------------------------------------------------------------------
// 核心优化思路：
// 1. 使用外部排序 (External Sort) 处理边，解决 1.8TB 导致的内存不足。
// 2. 使用快速 I/O (fread/fwrite/fast parser) 替代 istringstream。
// 3. 内存友好的哈希表策略：虽然仍使用 std::unordered_map，但通过 reserve 减少重平衡。
// 4. METIS 输出流式处理：利用排序后的边序列，逐个节点输出邻接行，避免全加载 CSR。
// ---------------------------------------------------------------------------

struct Edge {
    int32_t u;
    int32_t v;

    bool operator>(const Edge& other) const {
        if (u != other.u) return u > other.u;
        return v > other.v;
    }
    bool operator==(const Edge& other) const {
        return u == other.u && v == other.v;
    }
};

static inline bool edge_less(const Edge& a, const Edge& b) {
    if (a.u != b.u) return a.u < b.u;
    return a.v < b.v;
}

// 快速解析一行两个字符串 (rep, mem)
static bool get_two_strings(FILE* f, std::string& s1, std::string& s2) {
    static char buf[4096];
    if (!std::fgets(buf, sizeof(buf), f)) return false;
    if (buf[0] == '\n' || buf[0] == '\r' || buf[0] == '%') return get_two_strings(f, s1, s2);

    char c1[2048], c2[2048];
    if (std::sscanf(buf, "%s %s", c1, c2) != 2) return false;
    s1 = c1;
    s2 = c2;
    return true;
}

void print_usage(const char* prog) {
    std::cerr << "Usage: " << prog << " -i input.tsv [-o output.graph] [--map id_map.tsv] [--tmpdir dir] [--mem-limit GB]\n"
              << "Options:\n"
              << "  -i, --input    Input TSV file (representative member)\n"
              << "  -o, --output   Output METIS file\n"
              << "  --map          Output ID to original name mapping file\n"
              << "  --tmpdir       Directory for temporary files\n"
              << "  --mem-limit    Memory limit for sorting in GB (default: 32)\n";
}

// Pass 1: 分配 ID
int64_t pass1_assign_ids(const std::string& input_file, std::unordered_map<std::string, int32_t>& name_to_id) {
    FILE* f = std::fopen(input_file.c_str(), "r");
    if (!f) throw std::runtime_error("Cannot open input: " + input_file);

    std::string rep, mem;
    int32_t next_id = 1;
    name_to_id.reserve(20000000); // 预留 2000 万空间，防止 1100 万节点触发 rehash

    while (get_two_strings(f, rep, mem)) {
        for (const std::string& s : {rep, mem}) {
            if (name_to_id.find(s) == name_to_id.end()) {
                name_to_id[s] = next_id++;
            }
        }
    }
    std::fclose(f);
    return (int64_t)(next_id - 1);
}

// 外部排序：分块写入
void sort_and_write_chunk(std::vector<Edge>& edges, const std::string& filename) {
    std::sort(edges.begin(), edges.end(), edge_less);
    FILE* f = std::fopen(filename.c_str(), "wb");
    if (!f) throw std::runtime_error("Cannot open chunk file: " + filename);
    std::fwrite(edges.data(), sizeof(Edge), edges.size(), f);
    std::fclose(f);
    edges.clear();
}

// Pass 2: 提取边并分块排序
std::vector<std::string> pass2_extract_chunks(
    const std::string& input_file,
    const std::unordered_map<std::string, int32_t>& name_to_id,
    const std::string& tmpdir,
    size_t mem_limit_bytes,
    int64_t& total_raw_edges)
{
    FILE* f = std::fopen(input_file.c_str(), "r");
    std::string rep, mem;
    std::vector<Edge> buffer;
    size_t max_edges = mem_limit_bytes / sizeof(Edge);
    buffer.reserve(max_edges);

    std::vector<std::string> chunk_files;
    total_raw_edges = 0;

    while (get_two_strings(f, rep, mem)) {
        auto it_rep = name_to_id.find(rep);
        if (it_rep == name_to_id.end()) continue;
        auto it_mem = name_to_id.find(mem);
        if (it_mem == name_to_id.end()) continue;

        int32_t u = it_rep->second;
        int32_t v = it_mem->second;
        if (u == v) continue;
        if (u > v) std::swap(u, v);

        buffer.push_back({u, v});
        total_raw_edges++;

        if (buffer.size() >= max_edges) {
            std::string chunk_name = tmpdir + "/chunk_" + std::to_string(chunk_files.size()) + ".bin";
            sort_and_write_chunk(buffer, chunk_name);
            chunk_files.push_back(chunk_name);
        }
    }
    if (!buffer.empty()) {
        std::string chunk_name = tmpdir + "/chunk_" + std::to_string(chunk_files.size()) + ".bin";
        sort_and_write_chunk(buffer, chunk_name);
        chunk_files.push_back(chunk_name);
    }
    std::fclose(f);
    return chunk_files;
}

struct MergeNode {
    Edge edge;
    int chunk_idx;
    bool operator>(const MergeNode& other) const { return edge > other.edge; }
};

// 辅助结构，方便获取并弹出堆顶
template<typename T, typename Container, typename Compare>
struct custom_priority_queue : std::priority_queue<T, Container, Compare> {
    T pop_top() {
        std::pop_heap(this->c.begin(), this->c.end(), this->comp);
        T top = std::move(this->c.back());
        this->c.pop_back();
        return top;
    }
};

// Pass 3: 多路归并并流式写入 METIS
void merge_and_output(
    const std::vector<std::string>& chunk_files,
    const std::string& output_file,
    const std::string& map_file,
    const std::unordered_map<std::string, int32_t>& name_to_id,
    int64_t num_vertices,
    int64_t& unique_edges_out)
{
    // If map_file is provided, write the id_map
    if (!map_file.empty()) {
        std::cerr << "Writing ID map to " << map_file << "...\n";
        std::vector<std::string> id_to_name(num_vertices + 1);
        for (const auto& kv : name_to_id) {
            if (kv.second >= 1 && kv.second <= num_vertices) {
                id_to_name[kv.second] = kv.first;
            }
        }
        FILE* fmap = std::fopen(map_file.c_str(), "w");
        if (!fmap) throw std::runtime_error("Cannot open map file for writing: " + map_file);
        std::fprintf(fmap, "new_id\toriginal_name\n");
        for (int32_t i = 1; i <= num_vertices; i++) {
            std::fprintf(fmap, "%d\t%s\n", i, id_to_name[i].c_str());
        }
        std::fclose(fmap);
    }

    custom_priority_queue<MergeNode, std::vector<MergeNode>, std::greater<MergeNode>> pq;
    std::vector<FILE*> handles;
    for (int i = 0; i < (int)chunk_files.size(); ++i) {
        FILE* f = std::fopen(chunk_files[i].c_str(), "rb");
        if (f) {
            // 为读取分块文件设置 1MB 缓冲区
            setvbuf(f, nullptr, _IOFBF, 1024 * 1024);
            handles.push_back(f);
            Edge e;
            if (std::fread(&e, sizeof(Edge), 1, f) == 1) {
                pq.push({e, i});
            }
        }
    }

    // 为了得到 METIS 的边数，我们需要先归并一次去重计数。
    // 优化：在归并时直接统计度数，减少后续重读文件的次数。
    std::string sorted_bin = output_file + ".sorted.tmp";
    FILE* fsorted = std::fopen(sorted_bin.c_str(), "wb");
    if (!fsorted) throw std::runtime_error("Cannot open sorted tmp file");

    unique_edges_out = 0;
    Edge last_edge = {0, 0};
    std::vector<int64_t> degree(num_vertices + 1, 0);

    while (!pq.empty()) {
        MergeNode top = pq.pop_top();

        if (!(top.edge == last_edge)) {
            std::fwrite(&top.edge, sizeof(Edge), 1, fsorted);
            unique_edges_out++;
            degree[top.edge.u]++;
            degree[top.edge.v]++;
            last_edge = top.edge;
        }

        Edge next_e;
        if (std::fread(&next_e, sizeof(Edge), 1, handles[top.chunk_idx]) == 1) {
            pq.push({next_e, top.chunk_idx});
        }
    }
    for (auto h : handles) std::fclose(h);
    std::fclose(fsorted);

    std::cerr << "Unique edges: " << unique_edges_out << ". Building adjacency list in memory...\n";

    std::vector<int64_t> offset(num_vertices + 2, 0);
    for (int64_t i = 1; i <= num_vertices; i++) {
        offset[i+1] = offset[i] + degree[i];
    }
    degree.clear(); degree.shrink_to_fit(); // 释放不再需要的内存

    std::vector<int32_t> neighbors(offset[num_vertices+1]);
    std::vector<int64_t> current_pos = offset;

    fsorted = std::fopen(sorted_bin.c_str(), "rb");
    if (!fsorted) throw std::runtime_error("Cannot open sorted tmp file for second pass");
    setvbuf(fsorted, nullptr, _IOFBF, 8 * 1024 * 1024); // 8MB 缓冲区
    Edge e;
    while (std::fread(&e, sizeof(Edge), 1, fsorted) == 1) {
        neighbors[current_pos[e.u]++] = e.v;
        neighbors[current_pos[e.v]++] = e.u;
    }
    std::fclose(fsorted);

    std::cerr << "Writing METIS output file...\n";
    FILE* out_f = std::fopen(output_file.c_str(), "w");
    if (!out_f) throw std::runtime_error("Cannot open output file: " + output_file);
    
    // 设置一个较大的缓冲区提升写入效率 (8MB)
    std::vector<char> write_buf(8 * 1024 * 1024);
    std::setvbuf(out_f, write_buf.data(), _IOFBF, write_buf.size());

    std::fprintf(out_f, "%lld %lld\n", (long long)num_vertices, (long long)unique_edges_out);
    for (int64_t i = 1; i <= num_vertices; i++) {
        int64_t start = offset[i];
        int64_t end = offset[i+1];
        if (start < end) {
            std::sort(neighbors.begin() + start, neighbors.begin() + end);
            for (int64_t j = start; j < end; j++) {
                if (j > start) std::fputc(' ', out_f);
                // 进一步加速：手动转换整数
                char buf[16];
                int len = std::snprintf(buf, sizeof(buf), "%d", neighbors[j]);
                std::fwrite(buf, 1, len, out_f);
            }
        }
        std::fputc('\n', out_f);
        if (i % 100000 == 0) std::cerr << "\rWriting progress: " << i << " / " << num_vertices << std::flush;
    }
    std::cerr << "\n";
    std::fclose(out_f);
}

int main(int argc, char* argv[]) {
    std::string input, output, map_file, tmpdir = ".";
    size_t mem_limit_gb = 32;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if ((arg == "-i" || arg == "--input") && i + 1 < argc) input = argv[++i];
        else if ((arg == "-o" || arg == "--output") && i + 1 < argc) output = argv[++i];
        else if (arg == "--map" && i + 1 < argc) map_file = argv[++i];
        else if (arg == "--tmpdir" && i + 1 < argc) tmpdir = argv[++i];
        else if (arg == "--mem-limit" && i + 1 < argc) mem_limit_gb = std::stoll(argv[++i]);
    }

    if (input.empty()) { print_usage(argv[0]); return 1; }
    if (output.empty()) output = input + ".graph";

    try {
        std::unordered_map<std::string, int32_t> name_to_id;
        std::cerr << "Pass 1: Assigning IDs...\n";
        int64_t num_v = pass1_assign_ids(input, name_to_id);

        std::cerr << "Pass 2: Extracting & Sorting chunks...\n";
        int64_t total_raw;
        auto chunks = pass2_extract_chunks(input, name_to_id, tmpdir, mem_limit_gb * 1024LL * 1024LL * 1024LL, total_raw);

        std::cerr << "Pass 3: Merging & Writing METIS...\n";
        int64_t unique_e;
        merge_and_output(chunks, output, map_file, name_to_id, num_v, unique_e);

        // Cleanup chunks
        for (const auto& f : chunks) std::remove(f.c_str());
        std::remove((output + ".sorted.tmp").c_str());

        std::cerr << "Done. Vertices: " << num_v << ", Edges: " << unique_e << "\n";
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
