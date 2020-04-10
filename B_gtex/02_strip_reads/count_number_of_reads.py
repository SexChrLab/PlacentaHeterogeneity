# In this script, I want to count the number of reads after stripping the reads

import os
import sys
import gzip

id = sys.argv[1]
outfile = open(sys.argv[2], 'w')

id_rg = []

rg_fn = os.path.join('/scratch/tphung3/PlacentaSexDiff/B_gtex/02_strip_reads/strip_reads', id, 'fastq', id + '.full_rg.list')
with open(rg_fn, 'r') as f:
    for line in f:
        id_rg.append(line.rstrip('\n'))

for i in id_rg:
    for j in (1, 2):
        read_fn = os.path.join('/scratch/tphung3/PlacentaSexDiff/B_gtex/02_strip_reads/strip_reads', id, 'fastq', id + '_' + i + '_' + str(j) + '.fastq.gz')
        read_nreads = 0
        with gzip.open(read_fn, 'r') as f:
            for line in f:
                read_nreads += 1
        out = [id, i, 'read_' + str(j), str(read_nreads)]
        print ('\t'.join(out), file=outfile)
