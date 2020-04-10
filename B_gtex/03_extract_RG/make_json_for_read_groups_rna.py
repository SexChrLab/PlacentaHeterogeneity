# In this script, I will generate the information for each read group in a json format so I can use this as part of the config file for the snakemake pipeline

import os
import json
from collections import defaultdict

with open('/scratch/tphung3/PlacentaSexDiff/B_gtex/04_process_rna/process_rna_config.json') as f:
    data = json.load(f)

rna_read_group_identifier = defaultdict(list)
rna_samples = {"rna_samples": defaultdict(list)}

for filename in os.listdir('/scratch/tphung3/PlacentaSexDiff/B_gtex/03_extract_RG/read_groups/rna/'):
    with open(os.path.join('/scratch/tphung3/PlacentaSexDiff/B_gtex/03_extract_RG/read_groups/rna/', filename), 'r') as file:
        filename_id = filename.split('_')[0]
        for line in file:
            read_group_info = {}
            items = line.rstrip('\n').split('\t')
            # print (items)
            for i in items:
                if i.startswith('ID'):
                    id = i.split(':')[1]
                if i.startswith('SM'):
                    sm = i.split(':')[1]
                if i.startswith('LB'):
                    lb = i.split(':')[1]
                if i.startswith('PU'):
                    pu = i.split(':')[1]
                if i.startswith('PL'):
                    pl = i.split(':')[1]

            key = filename_id + '_' + id
            fq_path = '/scratch/tphung3/PlacentaSexDiff/B_gtex/02_strip_reads/strip_reads/' + filename_id + '/fastq/'
            fq_1 = filename_id + '_' + id + '_1.fastq.gz'
            fq_2 = filename_id + '_' + id + '_2.fastq.gz'

            rna_read_group_identifier["rna_read_group_identifier"].append(key)

            read_group_info[key] = {
                'fq_path': fq_path,
                'fq_1': fq_1,
                'fq_2': fq_2,
                'ID': id,
                'SM': sm,
                'LB': lb,
                'PU': pu,
                'PL': pl
            }

            data.update(read_group_info)

            with open('/scratch/tphung3/PlacentaSexDiff/B_gtex/04_process_rna/process_rna_config.json', 'w') as f:
                json.dump(data, f)

            rna_samples["rna_samples"][filename_id].append(key)

data.update(rna_read_group_identifier)
data.update(rna_samples)

with open('/scratch/tphung3/PlacentaSexDiff/B_gtex/04_process_rna/process_rna_config.json', 'w') as f:
    json.dump(data, f)
