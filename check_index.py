lines = open('/Users/zhixy/CLionProjects/mmseqs_ccc/test/filtered_db.index').readlines()
ids = [int(l.split()[0]) for l in lines]
offsets = [int(l.split()[1]) for l in lines]
sorted_by_offset = all(offsets[i] <= offsets[i+1] for i in range(len(offsets)-1))
sorted_by_id = all(ids[i] <= ids[i+1] for i in range(len(ids)-1))
print('index sorted by offset:', sorted_by_offset)
print('index sorted by query_id:', sorted_by_id)
print('sample offsets:', offsets[:10])
print('sample ids:', ids[:10])
