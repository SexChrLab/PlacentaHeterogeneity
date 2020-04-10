# In this script, I want to count the number of reads after stripping the reads

import os
import sys
import gzip
import os.path
from os import path

id = sys.argv[1]

if path.exists(sys.argv[2]):
    outfile = open(sys.argv[2], 'a')
else:
    outfile = open(sys.argv[2], 'w')

id_rg = []

rg_fn = os.path.join('/scratch/tphung3/PlacentaSexDiff/B_gtex/02_strip_reads/strip_reads', id, 'fastq', id + '.full_rg.list')
with open(rg_fn, 'r') as f:
    for line in f:
        id_rg.append(line.rstrip('\n'))

for i in id_rg:
    read_1_fn = os.path.join('/scratch/tphung3/PlacentaSexDiff/B_gtex/02_strip_reads/strip_reads', id, 'fastq', id + '_' + i + '_' + str(1) + '.fastq.gz')
    read_1_nreads = 0
    with gzip.open(read_1_fn, 'r') as f:
        for line in f:
            read_1_nreads += 1
    read_2_fn = os.path.join('/scratch/tphung3/PlacentaSexDiff/B_gtex/02_strip_reads/strip_reads', id, 'fastq', id + '_' + i + '_' + str(2) + '.fastq.gz')
    read_2_nreads = 0
    with gzip.open(read_2_fn, 'r') as f:
        for line in f:
            read_2_nreads += 1
    if read_1_nreads != read_2_nreads:
        out_1 = [id, i, 'read_' + str(1), str(read_1_nreads)]
        print ('\t'.join(out_1), file=outfile)
        out_2 = [id, i, 'read_' + str(2), str(read_2_nreads)]
        print ('\t'.join(out_2), file=outfile)
