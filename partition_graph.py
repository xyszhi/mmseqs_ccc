#!/usr/bin/env python3
"""
用于巨型网络拆分的KaMinPar Python脚本
用法: python partition_graph.py <input.metis> <num_blocks> [-o OUTPUT] [-e EPS] [-t THREADS]
"""

import kaminpar
import argparse
import time
import sys
import os
from collections import Counter

def main():
    total_start = time.time()

    # --- 1. 解析命令行参数 ---
    parser = argparse.ArgumentParser(description="使用KaMinPar对巨型网络进行分区")
    parser.add_argument("input_file", help="输入的METIS格式图文件路径")
    parser.add_argument("k", type=int, help="需要划分成的区块数量 (k)")
    parser.add_argument("-o", "--output", default=None,
                        help="输出分区结果的文件路径 (默认: 与输入文件同目录，命名为 <basename>_partition.txt)")
    parser.add_argument("-e", "--epsilon", type=float, default=0.03,
                        help="允许的最大不平衡度，例如0.03表示3%% (默认: 0.03)")
    parser.add_argument("-t", "--threads", type=int, default=0,
                        help="使用的线程数，默认为0 (自动检测所有可用核心)")

    args = parser.parse_args()

    # 检查输入文件是否存在
    if not os.path.exists(args.input_file):
        print(f"错误: 输入文件 '{args.input_file}' 未找到。")
        sys.exit(1)

    # 默认输出路径与输入文件同目录
    if args.output is None:
        base, _ = os.path.splitext(args.input_file)
        args.output = base + "_partition.txt"

    # --- 2. 准备工作: 设置线程数和上下文 ---
    num_threads = args.threads if args.threads > 0 else os.cpu_count()
    print(f"检测到可用CPU核心数: {os.cpu_count()}")
    print(f"将使用 {num_threads} 个线程。")

    # 获取默认配置（这是Python绑定的标准用法）
    ctx = kaminpar.default_context()

    # 初始化KaMinPar实例
    print("正在初始化KaMinPar实例...")
    instance = kaminpar.KaMinPar(num_threads=num_threads, ctx=ctx)

    # --- 3. 加载图 ---
    print(f"正在从 '{args.input_file}' 加载图...")
    start_time = time.time()
    try:
        # 加载METIS格式的图
        graph = kaminpar.load_graph(args.input_file, kaminpar.GraphFileFormat.METIS, compress=False)
        print(f"图加载成功。")
        print(f"  - 节点数 (Vertices): {graph.n():,}")
        print(f"  - 边数 (Edges): {graph.m():,}")
    except Exception as e:
        print(f"加载图时发生错误: {e}")
        print("请检查文件格式是否为有效的METIS格式。")
        sys.exit(1)
    load_time = time.time() - start_time
    print(f"图加载耗时: {load_time:.2f} 秒。")

    # --- 4. 执行分区 ---
    print(f"\n开始将图划分为 {args.k} 个区块，允许的最大不平衡度为 {args.epsilon*100:.1f}%...")
    print("这个过程可能需要一些时间，请耐心等待...")
    partition_start = time.time()

    try:
        # compute_partition 返回一个数组，partition[i] 是节点 i 所属的区块ID (0-based)
        partition = instance.compute_partition(graph, k=args.k, eps=args.epsilon)
    except Exception as e:
        print(f"\n分区计算过程中发生错误: {e}")
        print("这可能是因为输入文件太大或内存不足。")
        sys.exit(1)

    partition_time = time.time() - partition_start
    print(f"分区计算完成！核心算法耗时: {partition_time:.2f} 秒。")

    # --- 5. 分析结果: 计算边割 (Edge Cut) ---
    edge_cut = kaminpar.edge_cut(graph, partition)
    m = graph.m()
    ratio = edge_cut / m if m > 0 else 0.0
    print(f"\n分区结果分析:")
    print(f"  - 总边割 (Edge Cut): {edge_cut:,}")
    print(f"  - 总边割比 (Edge Cut Ratio): {ratio:.6f}")

    # 简单统计各区块大小，检查负载均衡情况
    block_sizes = Counter(partition)
    actual_k = len(block_sizes)
    print(f"  - 区块数量: {actual_k}")
    if actual_k < args.k:
        print(f"  ⚠ 警告: 实际区块数 ({actual_k}) 少于请求的 k={args.k}")
    print(f"  - 最大区块大小: {max(block_sizes.values()):,}")
    print(f"  - 最小区块大小: {min(block_sizes.values()):,}")
    print(f"  - 平均区块大小: {sum(block_sizes.values()) / actual_k:.2f}")

    # --- 6. 保存分区结果 ---
    print(f"\n正在将分区结果保存到 '{args.output}'...")
    try:
        with open(args.output, 'w') as f:
            f.write(f"# KaMinPar Partition Result\n")
            f.write(f"# Input Graph: {args.input_file}\n")
            f.write(f"# Number of Blocks (k): {args.k}\n")
            f.write(f"# Epsilon: {args.epsilon}\n")
            f.write(f"# Edge Cut: {edge_cut}\n")
            f.write(f"# Nodes: {graph.n()}\n")
            f.write("# Block assignment for each node (one per line, 0-based index):\n")
            f.writelines(f"{b}\n" for b in partition)
        print(f"分区结果已成功保存。")

        # 同时保存一个区块摘要文件
        base, ext = os.path.splitext(args.output)
        summary_file = base + '_summary' + (ext if ext else '.txt')
        with open(summary_file, 'w') as f:
            f.write("Block ID\tSize\n")
            for block_id, size in sorted(block_sizes.items()):
                f.write(f"{block_id}\t{size}\n")
        print(f"区块摘要已保存到: {summary_file}")

    except Exception as e:
        print(f"保存结果时出错: {e}")

    total_time = time.time() - total_start
    print(f"\n脚本执行完毕。总耗时: {total_time:.2f} 秒。")

if __name__ == "__main__":
    main()
