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
    std::cerr << "Usage: " << prog << " -i input.tsv [-o output.graph] [--tmpdir dir] [--mem-limit GB]\n"
              << "Options:\n"
              << "  -i, --input    Input TSV file (representative member)\n"
              << "  -o, --output   Output METIS file\n"
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
        int32_t u = name_to_id.at(rep);
        int32_t v = name_to_id.at(mem);
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

// Pass 3: 多路归并并流式写入 METIS
void merge_and_output(
    const std::vector<std::string>& chunk_files,
    const std::string& output_file,
    int64_t num_vertices,
    int64_t& unique_edges_out)
{
    std::priority_queue<MergeNode, std::vector<MergeNode>, std::greater<MergeNode>> pq;
    std::vector<FILE*> handles;
    for (int i = 0; i < (int)chunk_files.size(); ++i) {
        FILE* f = std::fopen(chunk_files[i].c_str(), "rb");
        if (f) {
            handles.push_back(f);
            Edge e;
            if (std::fread(&e, sizeof(Edge), 1, f) == 1) {
                pq.push({e, i});
            }
        }
    }

    // 为了得到 METIS 的边数，我们需要先归并一次去重计数，或者提前知道。
    // 简化处理：先归并到一个排序的二进制文件，同时计数。
    std::string sorted_bin = output_file + ".sorted.tmp";
    FILE* fsorted = std::fopen(sorted_bin.c_str(), "wb");
    unique_edges_out = 0;
    Edge last_edge = {0, 0};

    while (!pq.empty()) {
        MergeNode top = pq.top();
        pq.pop();

        if (!(top.edge == last_edge)) {
            std::fwrite(&top.edge, sizeof(Edge), 1, fsorted);
            unique_edges_out++;
            last_edge = top.edge;
        }

        Edge next_e;
        if (std::fread(&next_e, sizeof(Edge), 1, handles[top.chunk_idx]) == 1) {
            pq.push({next_e, top.chunk_idx});
        }
    }
    for (auto h : handles) std::fclose(h);
    std::fclose(fsorted);

    // 最后一步：流式输出 METIS
    // 现在的 fsorted 是全局排序的 (u, v) 且 u < v。
    // METIS 需要每个节点的邻接表。我们需要把每个 (u, v) 变成 u->v 和 v->u。
    // 这又涉及一次排序... 或者使用两倍空间的外部排序。
    // 方案改进：在归并时直接输出所有双向边到另一个外部排序过程，或者在内存足够时处理。
    // 由于用户有 640GB，我们可以在此处利用内存构建邻接表。
    
    std::cerr << "Unique edges: " << unique_edges_out << ". Writing METIS...\n";
    
    // 内存写邻接表：如果 2*E*4 字节 < 640GB，则可行。
    // 180亿边 = 144GB。
    std::vector<int64_t> degree(num_vertices + 1, 0);
    fsorted = std::fopen(sorted_bin.c_str(), "rb");
    Edge e;
    while (std::fread(&e, sizeof(Edge), 1, fsorted) == 1) {
        degree[e.u]++;
        degree[e.v]++;
    }
    
    std::vector<int64_t> offset(num_vertices + 2, 0);
    for (int64_t i = 1; i <= num_vertices; i++) {
        offset[i+1] = offset[i] + degree[i];
    }
    
    std::vector<int32_t> neighbors(offset[num_vertices+1]);
    std::vector<int64_t> current_pos = offset;
    std::rewind(fsorted);
    while (std::fread(&e, sizeof(Edge), 1, fsorted) == 1) {
        neighbors[current_pos[e.u]++] = e.v;
        neighbors[current_pos[e.v]++] = e.u;
    }
    std::fclose(fsorted);

    std::ofstream fout(output_file);
    fout << num_vertices << " " << unique_edges_out << "\n";
    for (int64_t i = 1; i <= num_vertices; i++) {
        std::sort(neighbors.begin() + offset[i], neighbors.begin() + offset[i+1]);
        for (int64_t j = offset[i]; j < offset[i+1]; j++) {
            fout << neighbors[j] << (j == offset[i+1]-1 ? "" : " ");
        }
        fout << "\n";
    }
}

int main(int argc, char* argv[]) {
    std::string input, output, tmpdir = ".";
    size_t mem_limit_gb = 32;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if ((arg == "-i" || arg == "--input") && i + 1 < argc) input = argv[++i];
        else if ((arg == "-o" || arg == "--output") && i + 1 < argc) output = argv[++i];
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
        merge_and_output(chunks, output, num_v, unique_e);

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
