import os

db_prefix = '/Users/zhixy/CLionProjects/mmseqs_ccc/test/filtered_db'
num_splits = 56

sizes = [os.path.getsize(f'{db_prefix}.{i}') for i in range(num_splits)]
cum = [0] * (num_splits + 1)
for i in range(num_splits):
    cum[i+1] = cum[i] + sizes[i]

lines = open(f'{db_prefix}.index').readlines()

# Count seekg operations per file if we group by file
from collections import defaultdict
file_queries = defaultdict(list)
for l in lines:
    parts = l.split()
    qid, offset, length = int(parts[0]), int(parts[1]), int(parts[2])
    fid = 0
    for i in range(num_splits):
        if cum[i] <= offset < cum[i+1]:
            fid = i
            break
    file_queries[fid].append((offset - cum[fid], length, qid))

print('Queries per split file:')
for fid in sorted(file_queries.keys()):
    entries = file_queries[fid]
    offsets = sorted(e[0] for e in entries)
    # Check if sorted offsets are monotonically increasing (sequential reads possible)
    is_sorted = all(offsets[i] <= offsets[i+1] for i in range(len(offsets)-1))
    print(f'  file {fid}: {len(entries)} queries, local offsets sorted={is_sorted}')
