# In this script, we are generating the config json file for processing the exomes of placenta project
# - 12 female placentas from batch 1
# - 18 female placentas from batch 2
# - 18 female deciduas from batch 2

import json
from collections import defaultdict
import os
import gzip

data = {}

data['all_dna_samples'] = []
dna_samples_fastq_path = defaultdict(list)

with open('/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/female_sample_ids.csv', 'r') as f:
    for line in f:
        items = line.rstrip('\n').split(',')
        data['all_dna_samples'].append(items[0])
        dna_samples_fastq_path[items[0]].append(items[1])
        dna_samples_fastq_path[items[0]].append(items[2])

for sample in data['all_dna_samples']:
    read_group_info = {}
    if not '-' in sample:
        fq_path = '/data/mwilsons/lab_generated/placenta/fastqs/placentas/exome/'
        fq_1 = os.path.join(sample, 'Unaligned', dna_samples_fastq_path[sample][0])
        fq_2 = os.path.join(sample, 'Unaligned', dna_samples_fastq_path[sample][1])

        # find pu
        with gzip.open(os.path.join(fq_path, fq_1), 'rt', encoding="utf8", errors='ignore') as f:
            first_line = f.readline() #@K00175:45:HJC3YBBXX:6:1101:1742:1314 1:N:0:NGTCAA+NACAAG
            items = first_line.split(':')
            pu = items[2] + '.' + items[3]

        read_group_info[sample] = {
            'fq_path': fq_path,
            'fq_1': fq_1,
            'fq_2': fq_2,
            'ID': sample,
            'SM': sample,
            'LB': sample,
            'PU': pu,
            'PL': 'Illumina'
        }
    else:
        fq_path = '/data/mwilsons/lab_generated/placenta/fastqs/decidua_placenta/exome/'
        fq_1 = dna_samples_fastq_path[sample][0]
        fq_2 = dna_samples_fastq_path[sample][1]

        # find pu
        with gzip.open(os.path.join(fq_path, fq_1), 'rt', encoding="utf8", errors='ignore') as f:
            first_line = f.readline()
            items = first_line.split(':')
            pu = items[2] + '.' + items[3]

        read_group_info[sample] = {
            'fq_path': fq_path,
            'fq_1': fq_1,
            'fq_2': fq_2,
            'ID': sample,
            'SM': sample,
            'LB': sample,
            'PU': pu,
            'PL': 'Illumina'
        }

    data.update(read_group_info)

with open(
    '/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/process_dna_females_config.json',
    'w') as outfile:
    json.dump(data, outfile)