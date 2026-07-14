# mmseqs_ccc

基于 MMseqs2 全对全比对结果，对蛋白质序列进行连通分量聚类的 C++ 工具集，并提供图分区、子图拆分与可视化等后处理功能。

## 工具一览

| 工具 | 类型 | 功能简介 |
|------|------|---------|
| `mmseqs_ccc` | C++ 可执行程序 | 并行连通分量聚类，输出簇成员关系表 |
| `extract_cluster` | C++ 可执行程序 | 按簇代表序列名提取成员的比对结果 |
| `extract_cluster_to_metis2` | C++ 可执行程序 | 将簇比对结果转换为 METIS 图格式（基础版） |
| `extract_cluster_to_metis3` | C++ 可执行程序 | 将簇比对结果转换为 METIS 图格式（外部排序优化版，适合超大簇） |
| `check_metis_connectivity` | C++ 可执行程序 | 检查 METIS 图文件的连通性，统计连通分量 |
| `check_tsv_connectivity` | C++ 可执行程序 | 检查聚类 TSV 文件的连通性，统计连通分量 |
| `split_by_partition` | C++ 可执行程序 | 按 KaMinPar 分区结果将比对 TSV 拆分为各分区子文件 |
| `split_by_connectivity` | C++ 可执行程序 | 对比对 TSV 做连通分量分析，将各连通分量拆分为独立子文件 |
| `partition_graph.py` | Python 脚本 | 调用 KaMinPar 对 METIS 图进行图分区 |
| `sankey_kaminpar.py` | Python 脚本 | 生成 KaMinPar 多轮分区结果的 Sankey 可视化图 |

---

## 背景

本工具集用于处理 MMseqs2 `filterdb` 过滤后的比对结果数据库（`filtered_db`）：

1. `mmseqs_ccc` 通过并查集（Union-Find）算法将具有同源关系的蛋白质序列归并为连通分量，输出聚类结果（`filtered_db_clusters.tsv`）。
2. `extract_cluster` 读取聚类结果，按指定的簇代表序列名提取所有成员的比对记录，输出格式类似 `mmseqs convertalis --format-output`，支持 `qcov`、`tcov`、`bits` 等字段。
3. `extract_cluster_to_metis2` / `extract_cluster_to_metis3` 将指定簇的比对结果转换为 METIS 图格式，供 KaMinPar 等图分区工具使用。
4. `partition_graph.py` 调用 KaMinPar 对 METIS 图进行 k 路分区，输出分区结果文件。
5. `split_by_partition` 按分区结果将比对 TSV 拆分为各分区子文件，同时输出跨分区边记录。
6. `split_by_connectivity` 对比对 TSV 做连通分量分析，将各连通分量拆分为独立子文件，适合对分区后子图做进一步细分。

---

## 编译

```bash
cmake -DCMAKE_BUILD_TYPE=Release -B build_release
cmake --build build_release -j 8
```

编译后生成（均在 `build_release/` 目录下）：

- `mmseqs_ccc`
- `extract_cluster`
- `extract_cluster_to_metis2`
- `extract_cluster_to_metis3`
- `check_metis_connectivity`
- `check_tsv_connectivity`
- `split_by_partition`
- `split_by_connectivity`

---

## 典型工作流程

```
mmseqs 建库 → 全对全搜索 → 过滤
       ↓
  mmseqs_ccc（连通分量聚类）
       ↓
  extract_cluster_to_metis3（提取超大簇为 METIS 图）
       ↓
  partition_graph.py（KaMinPar 图分区）
       ↓
  split_by_partition（按分区拆分比对 TSV）
       ↓
  split_by_connectivity（对各分区子图再做连通分量拆分）
```

---

## C++ 程序详细说明

### mmseqs_ccc — 连通分量聚类

#### 运行

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

#### 输出

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

### extract_cluster — 提取簇比对结果

#### 运行

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
| `alnlen` | 比对长度 |
| `mismatch` | 错配数 |
| `gapopen` | gap 开启次数 |
| `qstart` | 查询起始位置（0-based） |
| `qend` | 查询终止位置（0-based） |
| `qlen` | 查询序列长度 |
| `qcov` | 查询覆盖度 |
| `tstart` | 目标起始位置（0-based） |
| `tend` | 目标终止位置（0-based） |
| `tlen` | 目标序列长度 |
| `tcov` | 目标覆盖度 |
| `evalue` | E 值 |
| `bits` | Bit score |

> **注意：** `alnlen`、`mismatch`、`gapopen` 在无 backtrace 的数据库（`mmseqs search` 未加 `-a`）中为近似值。若需精确值，请在 `mmseqs search` 时加 `-a` 参数。

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

### extract_cluster_to_metis2 / extract_cluster_to_metis3 — 提取簇为 METIS 图

将指定簇的比对结果转换为 METIS 图格式（`.metis` 文件），供 KaMinPar 等图分区工具使用。同时输出节点名称到整数 ID 的映射文件（`.metis.id_map`）。

- **`extract_cluster_to_metis2`**：基础版，适合中等规模的簇。
- **`extract_cluster_to_metis3`**：外部排序优化版，使用 `fread`/`fwrite` 快速 I/O 和外部归并排序，适合超大簇（数百万节点、TB 级比对数据）。

#### 运行

```bash
./build_release/extract_cluster_to_metis3 [OPTIONS]
```

**主要参数：**

| 参数 | 说明 |
|------|------|
| `--tsv <path>` | 输入比对 TSV 文件（`extract_cluster` 输出） |
| `--out <path>` | 输出 METIS 图文件路径（同时生成 `.id_map`） |

**输出文件：**

- `<out>.metis`：METIS 格式图文件，节点为序列，边为比对关系。
- `<out>.metis.id_map`：节点整数 ID 到序列名称的映射（制表符分隔）。

---

### check_metis_connectivity — 检查 METIS 图连通性

读取 METIS `.graph` 文件，使用并查集统计连通分量数量及各分量大小。

#### 运行

```bash
./build_release/check_metis_connectivity --input <graph.metis> [--verbose]
```

**参数说明：**

| 参数 | 说明 |
|------|------|
| `--input <path>` / `-i <path>` | 输入 METIS 图文件（必填） |
| `--verbose` | 输出所有连通分量的大小分布（降序排列） |

**输出：**顶点数、边数、连通分量数、最大分量大小；`--verbose` 时输出完整分量大小分布。

---

### check_tsv_connectivity — 检查 TSV 聚类连通性

直接读取聚类 TSV 文件（`rep<TAB>member` 格式），使用并查集统计连通分量，无需先转换为 METIS 格式。

#### 运行

```bash
./build_release/check_tsv_connectivity --input <clusters.tsv> [--clique] [--verbose]
```

**参数说明：**

| 参数 | 说明 |
|------|------|
| `--input <path>` / `-i <path>` | 输入聚类 TSV 文件（必填） |
| `--clique` | 使用全连接拓扑（默认为星形拓扑：代表序列与每个成员相连） |
| `--verbose` | 输出连通分量大小分布（前 20 个） |

**退出码：** `0` = 单一连通分量，`2` = 多个连通分量，`1` = 错误。

---

### split_by_partition — 按分区拆分比对 TSV

读取 KaMinPar 分区结果（`partition.txt`）和节点 ID 映射（`.metis.id_map`），将比对 TSV 按分区拆分为各分区子文件。查询和目标属于不同分区的行输出到跨分区文件（`_cross.tsv`）。

#### 运行

```bash
./build_release/split_by_partition \
    --tsv <input.tsv> \
    --idmap <graph.metis.id_map> \
    --partition <partition.txt> \
    [--out-dir <dir>]
```

**参数说明：**

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--tsv <path>` | 输入比对 TSV 文件（必填） |
| `--idmap <path>` | METIS id_map 文件（必填） |
| `--partition <path>` | KaMinPar 分区结果文件（必填） |
| `--out-dir <path>` | 输出目录 | 与 `--tsv` 同目录 |

**输出文件：**

- `<tsv_basename>_K0.tsv`、`<tsv_basename>_K1.tsv`、…：各分区内部比对记录。
- `<tsv_basename>_cross.tsv`：跨分区比对记录。

---

### split_by_connectivity — 按连通分量拆分比对 TSV

对比对 TSV 文件做连通分量分析（两遍扫描），将每个连通分量的比对记录拆分为独立子文件。支持多线程并行处理。

#### 运行

```bash
./build_release/split_by_connectivity \
    --input <input.tsv> \
    --output <out_dir> \
    [--sep <sep>] \
    [--threads <n>] \
    [--prefix <prefix>]
```

**参数说明：**

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `-i` / `--input <path>` | 输入比对 TSV 文件（必填） |
| `-o` / `--output <dir>` | 输出目录（必填） |
| `-s` / `--sep <char>` | 字段分隔符 | `\t` |
| `-t` / `--threads <n>` | 线程数 | 硬件并发数 |
| `-p` / `--prefix <str>` | 输出文件名前缀 | `component_` |

**输出：** 每个连通分量生成一个独立的 TSV 文件，文件名为 `<prefix><component_id>.tsv`。

---

## Python 脚本详细说明

### partition_graph.py — KaMinPar 图分区

调用 Python KaMinPar 库对 METIS 格式图进行 k 路分区，支持多个 k 值批量运行，输出分区结果文件和分区摘要。

#### 运行

```bash
python partition_graph.py <input.metis> <num_blocks> [-o OUTPUT] [-e EPS] [-t THREADS]
```

**参数说明：**

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `input.metis` | 输入 METIS 图文件（必填） |
| `num_blocks` | 分区数 k（必填） |
| `-o` / `--output <path>` | 输出分区结果文件路径 | `<input>.partition_k<k>.txt` |
| `-e` / `--epsilon <float>` | 负载均衡容忍度（imbalance factor） | `0.03` |
| `-t` / `--threads <n>` | 线程数 | 全部可用核心 |

**输出文件：**

- `partition.txt`：每行一个整数，表示对应节点所属的分区编号（0-based）。
- `partition_summary.txt`：分区摘要，包含各分区节点数统计。

**依赖：** `kaminpar`（Python 包）

---

### sankey_kaminpar.py — KaMinPar 分区 Sankey 可视化

读取多轮 KaMinPar 分区结果（不同 k 值），生成展示节点在各轮分区间流向的 Sankey 图，输出为 PNG 图片。

#### 运行

```bash
python sankey_kaminpar.py
```

> 脚本内部硬编码了输入路径和参数（k 值范围、阈值等），使用前需根据实际路径修改脚本中的配置。

**依赖：** `plotly`、`kaleido`

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

# 4. 连通分量聚类
./build_release/mmseqs_ccc /data/filtered_db /data/protein_db.lookup 56

# 5. 提取超大簇为 METIS 图
./build_release/extract_cluster_to_metis3 \
    --tsv cluster_WP001.tsv \
    --out cluster_WP001_graph.metis

# 6. KaMinPar 图分区
python partition_graph.py cluster_WP001_graph.metis 10

# 7. 按分区拆分比对 TSV
./build_release/split_by_partition \
    --tsv cluster_WP001.tsv \
    --idmap cluster_WP001_graph.metis.id_map \
    --partition partition.txt

# 8. 对各分区子图做连通分量拆分
./build_release/split_by_connectivity \
    --input cluster_WP001_K0.tsv \
    --output ./components/
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
- Python 3（`partition_graph.py`、`sankey_kaminpar.py`）
- Python 包：`kaminpar`、`plotly`、`kaleido`
