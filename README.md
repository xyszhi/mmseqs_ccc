# mmseqs_ccc

基于 MMseqs2 全对全比对结果，对蛋白质序列进行连通分量聚类的 C++ 程序。

## 背景

本程序用于处理 MMseqs2 `filterdb` 过滤后的比对结果数据库（`filtered_db`），通过并查集（Union-Find）算法将具有同源关系的蛋白质序列归并为连通分量，输出聚类结果。

## 使用方法

### 编译

```bash
cmake -DCMAKE_BUILD_TYPE=Release -B build_release
cmake --build build_release --target mmseqs_ccc -j 8
```

### 运行

```bash
./build_release/mmseqs_ccc <filtered_db前缀> <protein_db.lookup路径> <线程数>
```

**参数说明：**

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `filtered_db前缀` | filterdb 输出文件的路径前缀（不含 `.0`、`.index` 等后缀） | `test/filtered_db` |
| `protein_db.lookup路径` | MMseqs2 lookup 文件，提供序列 ID 到序列名称的映射 | `test/protein_db.lookup` |
| `线程数` | 并行处理的线程数，建议设为服务器核心数 | `8` |

**示例：**

```bash
# 后台运行，输出日志
nohup ./build_release/mmseqs_ccc \
    /data/filtered_db \
    /data/protein_db.lookup \
    56 \
    > clustering.log 2>&1 &

tail -f clustering.log
```

## 输入文件格式

### filtered_db（MMseqs2 分片数据库）

由 MMseqs2 `filterdb` 命令生成，包含以下文件：

- `filtered_db.0`、`filtered_db.1`、…、`filtered_db.N`：数据分片，每个分片为二进制格式，内部按查询序列分块存储比对结果。
- `filtered_db.index`：全局索引文件，每行三列（制表符分隔）：
  ```
  query_id    global_offset    block_length
  ```

每个数据块内，每行为一条比对记录（制表符分隔）：

```
target_id  seq_id  aln_len  mismatches  gap_opens  q_start  q_end  t_start  t_end  evalue  bitscore
```

### protein_db.lookup

MMseqs2 lookup 文件，每行两列（制表符分隔）：

```
sequence_id    sequence_name
```

用于将内部整数 ID 映射回原始序列名称（如 `WP_000001234.1`）。

## 算法设计

### 总体流程

```
读取 lookup 文件
       ↓
统计分片文件数量，计算各分片累计偏移量
       ↓
读取 index 文件，将每条查询分配到对应分片
       ↓
多线程并行读取各分片，提取比对边
       ↓
加锁批量写入并查集（Union-Find）
       ↓
收集连通分量，统计并输出聚类结果
```

### 1. 分片文件映射

MMseqs2 的 `filterdb` 输出使用全局偏移量（`global_offset`）索引所有查询，但数据实际分散在多个分片文件中。程序通过以下方式将全局偏移量映射到具体分片：

1. 依次用 POSIX `stat()` 检测 `filtered_db.0`、`filtered_db.1`、… 是否存在，统计分片总数 `num_splits`。
2. 记录各分片文件大小，计算累计偏移量数组 `cum_sizes`：
   ```
   cum_sizes[0] = 0
   cum_sizes[i+1] = cum_sizes[i] + size(filtered_db.i)
   ```
3. 对 index 中每条记录的 `global_offset`，用二分查找（`std::upper_bound`）定位所属分片 `fid`，并计算局部偏移量：
   ```
   local_offset = global_offset - cum_sizes[fid]
   ```

### 2. 顺序 I/O 优化

将每个分片内的查询按 `local_offset` 升序排序后再读取，确保对每个分片文件的访问是严格顺序的，最大化磁盘 I/O 效率（对 HDD 尤为重要）。

### 3. 并行处理

使用 C++11 `std::thread` 实现工作线程池：

- 主线程创建 `num_threads` 个工作线程。
- 各线程通过原子计数器 `next_file`（`std::atomic<int>`）竞争获取下一个待处理的分片编号，实现动态负载均衡。
- 每个线程独立打开并顺序读取其负责的分片文件，将解析出的 `(query_id, target_id)` 边对收集到线程本地的 `edges` 向量中。
- 处理完一个分片后，线程持有 `std::mutex` 锁，将本地 `edges` 批量写入全局并查集，然后释放锁，继续处理下一个分片。

这种"先本地收集、再批量加锁"的策略显著减少了锁竞争次数。

### 4. 并查集（Union-Find / DSU）

使用路径压缩（路径减半）+ 按秩合并的并查集，时间复杂度接近 O(α(n))（反阿克曼函数，实际近似 O(1)）：

```cpp
struct DSU {
    std::vector<int> parent, rank_;

    int find(int x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]]; // 路径减半
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
```

并查集大小为 `id2name.size()`（即 lookup 文件中的序列总数），初始化时每个节点的父节点指向自身。

### 5. 连通分量收集与输出

所有分片处理完毕后：

1. 遍历 lookup 文件中所有已知序列 ID，调用 `dsu.find(id)` 获取其根节点，按根节点分组，得到所有连通分量。
2. 统计聚类数量、最大聚类大小、单例数量。
3. 将结果写入 `<filtered_db前缀>_clusters.tsv`，格式为制表符分隔的两列：
   ```
   cluster_representative_name    member_name
   ```
   每个序列恰好出现在一行的第二列，第一列为其所在聚类的代表序列（根节点对应的序列名称）。

## 输出文件

`filtered_db_clusters.tsv`（制表符分隔）：

```
WP_000001234.1    WP_000001234.1
WP_000001234.1    WP_000005678.1
WP_000001234.1    WP_000009012.1
WP_000099999.1    WP_000099999.1
...
```

- 每个蛋白质恰好出现一次（在第二列）。
- 同一聚类的所有成员共享相同的第一列（代表序列）。
- 单例蛋白质的两列相同。

## 资源需求

以 17,710,805 个蛋白、56 个切片（每片 ~38 GB，共 ~2.1 TB）为例：

| 资源 | 预估用量 |
|------|---------|
| 内存 | ~10–30 GB（峰值） |
| 运行时间（NVMe/并行存储） | ~10–30 分钟 |
| 运行时间（HDD） | ~1–2 小时 |
| 输出文件大小 | ~900 MB |

## 依赖

- C++14 或更高标准
- POSIX 系统（Linux/macOS）
- CMake 3.10+
- pthreads（通过 `find_package(Threads)` 链接）

## 上游数据生成流程

本程序处理的数据由以下 MMseqs2 流程生成：

```bash
# 1. 建库
mmseqs createdb bacterial.fasta protein_db --threads 8

# 2. 全对全搜索
mmseqs search protein_db protein_db result_db tmp \
    -s 7.0 --min-seq-id 0.3 -c 0.5 --min-aln-len 60 \
    --cov-mode 0 -e 0.001 --max-seqs 10000 --threads 8

# 3. 过滤（保留 E值 ≤ 1e-10 的比对）
mmseqs filterdb result_db filtered_db \
    --comparison-operator le --filter-column 4 \
    --comparison-value 1e-10 --threads 8
```
