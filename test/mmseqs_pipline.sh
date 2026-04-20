#!/usr/bin/env sh

mmseqs createdb --threads 8 bacterial.fasta protein_db
mmseqs search protein_db protein_db result_db tmp -s 7.0 --min-seq-id 0.3 -c 0.5 --min-aln-len 60 --cov-mode 0 --threads 8 -e 0.001 --max-seqs 10000
mmseqs filterdb result_db filtered_db --comparison-operator le --filter-column 4 --comparison-value 1e-10 --threads 8