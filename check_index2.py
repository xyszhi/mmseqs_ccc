import os

db_prefix = '/Users/zhixy/CLionProjects/mmseqs_ccc/test/filtered_db'
num_splits = 56

# Get cumulative sizes
sizes = [os.path.getsize(f'{db_prefix}.{i}') for i in range(num_splits)]
cum = [0] * (num_splits + 1)
for i in range(num_splits):
    cum[i+1] = cum[i] + sizes[i]

lines = open(f'{db_prefix}.index').readlines()

# For each query, determine which file it falls in
file_ids = []
for l in lines[:200]:
    parts = l.split()
    qid, offset, length = int(parts[0]), int(parts[1]), int(parts[2])
    fid = 0
    for i in range(num_splits):
        if cum[i] <= offset < cum[i+1]:
            fid = i
            break
    file_ids.append(fid)

print('file_id for first 20 queries:', file_ids[:20])
# Check pattern
print('unique file_ids in first 200:', sorted(set(file_ids)))
# Check if it's round-robin
print('pattern (file_id % 56 == qid % 56)?', all(file_ids[i] == i % num_splits for i in range(len(file_ids))))
