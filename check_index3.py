import os

db_prefix = '/Users/zhixy/CLionProjects/mmseqs_ccc/test/filtered_db'
num_splits = 56

sizes = [os.path.getsize(f'{db_prefix}.{i}') for i in range(num_splits)]
cum = [0] * (num_splits + 1)
for i in range(num_splits):
    cum[i+1] = cum[i] + sizes[i]

lines = open(f'{db_prefix}.index').readlines()

# Find where file_id changes
prev_fid = 0
transitions = []
for idx, l in enumerate(lines):
    parts = l.split()
    qid, offset = int(parts[0]), int(parts[1])
    fid = 0
    for i in range(num_splits):
        if cum[i] <= offset < cum[i+1]:
            fid = i
            break
    if fid != prev_fid:
        transitions.append((idx, qid, fid))
        prev_fid = fid

print('File transitions (index_line, query_id, file_id):')
for t in transitions[:20]:
    print(t)
print('Total transitions:', len(transitions))
