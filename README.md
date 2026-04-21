# mmseqs_ccc

基于 MMseqs2 全对全比对结果，对蛋白质序列进行连通分量聚类的 C++ 工具集，包含两个独立程序：

- **`mmseqs_ccc`**：并行连通分量聚类，输出簇成员关系表
- **`extract_cluster`**：按簇代表序列名提取成员的比对结果，支持自定义输出字段

## 背景

本工具集用于处理 MMseqs2 `filterdb` 过滤后的比对结果数据库（`filtered_db`）：

1. `mmseqs_ccc` 通过并查集（Union-Find）算法将具有同源关系的蛋白质序列归并为连通分量，输出聚类结果（`filtered_db_clusters.tsv`）。
2. `extract_cluster` 读取聚类结果，按指定的簇代表序列名提取所有成员的比对记录，输出格式类似 `mmseqs convertalis --format-output`，支持 `qcov`、`tcov`、`bits` 等字段。

## 编译

```bash
cmake -DCMAKE_BUILD_TYPE=Release -B build_release
cmake --build build_release -j 8
```

编译后生成：
- `build_release/mmseqs_ccc`
- `build_release/extract_cluster`

---

## mmseqs_ccc — 连通分量聚类

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
nohup ./build_release/mmseqs_ccc \
    /data/filtered_db \
    /data/protein_db.lookup \
    56 \
    > clustering.log 2>&1 &

tail -f clustering.log
```

### 输出

`<filtered_db前缀>_clusters.tsv`（制表符分隔，两列）：

```
cluster_rep    member
WP_000001234.1    WP_000001234.1
WP_000001234.1    WP_000005678.1
WP_000001234.1    WP_000009012.1
WP_000099999.1    WP_000099999.1
...
```

- 每个蛋白质恰好出现一次（在第二列）。
- 同一聚类的所有成员共享相同的第一列（代表序列）。
- 单例蛋白质的两列相同。

---

## extract_cluster — 提取簇比对结果

### 运行

```bash
./build_release/extract_cluster [OPTIONS] <rep_name> [rep_name2 ...]
```

**参数说明：**

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `<rep_name>` | 簇代表序列名（`_clusters.tsv` 第一列），可指定多个 | 必填 |
| `--db <path>` | filtered_db 路径前缀 | `test/filtered_db` |
| `--lookup <path>` | MMseqs2 lookup 文件路径 | `<db目录>/protein_db.lookup` |
| `--clusters <path>` | 聚类结果 TSV 文件路径 | `<db前缀>_clusters.tsv` |
| `--format <fmt>` | 输出字段，逗号分隔（见下表） | 见下方默认值 |
| `--sep <char>` | 输出字段分隔符 | `\t`（制表符） |
| `--out <path>` | 输出文件路径 | 标准输出 |
| `--header` | 输出表头行 | 否 |

**默认 `--format`：**

```
query,target,fident,alnlen,mismatch,gapopen,qstart,qend,qlen,qcov,tstart,tend,tlen,tcov,evalue,bits
```

**可用字段：**

| 字段 | 说明 |
|------|------|
| `query` | 查询序列名 |
| `target` | 目标序列名 |
| `fident` | 序列一致性（0–1） |
| `alnlen` | 比对长度（`max(|qend-qstart|, |tend-tstart|) + 1`，有 backtrace 时用 cigar 精确计算） |
| `mismatch` | 错配数（无 backtrace 时为近似值） |
| `gapopen` | gap 开启次数（无 backtrace 时为 0） |
| `qstart` | 查询起始位置（0-based） |
| `qend` | 查询终止位置（0-based） |
| `qlen` | 查询序列长度 |
| `qcov` | 查询覆盖度（`(qend-qstart)/qlen`） |
| `tstart` | 目标起始位置（0-based） |
| `tend` | 目标终止位置（0-based） |
| `tlen` | 目标序列长度 |
| `tcov` | 目标覆盖度（`(tend-tstart)/tlen`） |
| `evalue` | E 值 |
| `bits` | Bit score |

> **注意：** `alnlen`、`mismatch`、`gapopen` 在无 backtrace 的数据库（`mmseqs search` 未加 `-a`）中为近似值：
> - `alnlen` = `max(|qend-qstart|, |tend-tstart|) + 1`（与 mmseqs2 `Matcher::computeAlnLength()` 一致）
> - `mismatch` = `min(|qend-qstart|, |tend-tstart|) × (1 - fident)`（与 mmseqs2 `convertalignments.cpp` 一致）
> - `gapopen` = 0
>
> 若需精确值，请在 `mmseqs search` 时加 `-a` 参数生成带 backtrace 的数据库。

**示例：**

```bash
# 提取单个簇，输出到文件，带表头
./build_release/extract_cluster \
    --db /data/filtered_db \
    --lookup /data/protein_db.lookup \
    --clusters /data/filtered_db_clusters.tsv \
    --header \
    --out cluster_WP001.tsv \
    "WP_000001234.1"

# 同时提取多个簇，自定义字段
./build_release/extract_cluster \
    --db /data/filtered_db \
    --format "query,target,fident,evalue,qcov,tcov,bits" \
    "WP_000001234.1" "WP_000099999.1"
```

---

## 上游数据生成流程

```bash
# 1. 建库
mmseqs createdb bacterial.fasta protein_db

# 2. 全对全搜索（加 -a 可获得精确 alnlen/mismatch/gapopen）
mmseqs search protein_db protein_db result_db tmp \
    -s 7.0 --min-seq-id 0.3 -c 0.5 --min-aln-len 60 \
    --cov-mode 0 -e 0.001 --max-seqs 10000 --threads 56 -a

# 3. 过滤（保留 E 值 ≤ 1e-10 的比对）
mmseqs filterdb result_db filtered_db \
    --comparison-operator le --filter-column 4 \
    --comparison-value 1e-10 --threads 56

# 4. 聚类
./build_release/mmseqs_ccc /data/filtered_db /data/protein_db.lookup 56

# 5. 提取指定簇的比对结果
./build_release/extract_cluster \
    --db /data/filtered_db \
    --header \
    "WP_000001234.1"
```

---

## 输入文件格式

### filtered_db（MMseqs2 分片数据库）

由 MMseqs2 `filterdb` 命令生成，包含以下文件：

- `filtered_db.0`、`filtered_db.1`、…、`filtered_db.N`：数据分片，二进制格式，内部按查询序列分块存储比对结果。
- `filtered_db.index`：全局索引文件，每行三列（制表符分隔）：
  ```
  query_id    global_offset    block_length
  ```
- `filtered_db.dbtype`：数据库类型标识文件。

每个数据块内，每行为一条比对记录（制表符分隔，10 或 11 列）：

```
target_id  score  fident  evalue  qstart  qend  qlen  tstart  tend  tlen  [cigar]
```

其中 `cigar` 列仅在 `mmseqs search -a` 时存在。

### protein_db.lookup

MMseqs2 lookup 文件，每行两列（制表符分隔）：

```
sequence_id    sequence_name
```

用于将内部整数 ID 映射回原始序列名称。

---

## 算法设计（mmseqs_ccc）

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
无锁并查集（atomic CAS）合并连通分量
       ↓
收集连通分量，统计并输出聚类结果
```

### 分片文件映射

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

### 顺序 I/O 优化

将每个分片内的查询按 `local_offset` 升序排序后再读取，确保对每个分片文件的访问是严格顺序的，最大化磁盘 I/O 效率（对 HDD 尤为重要）。

### 并行处理

使用 C++11 `std::thread` 实现工作线程池：

- 各线程通过原子计数器 `next_file`（`std::atomic<int>`）竞争获取下一个待处理的分片编号，实现动态负载均衡。
- 并查集使用无锁原子 CAS 操作（`compare_exchange`）实现并发合并，无需全局锁。

### 并查集（Union-Find / DSU）

使用路径压缩（路径减半）+ 按秩合并的无锁并查集（基于 `std::atomic<int>`），时间复杂度接近 O(α(n))。

---

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
