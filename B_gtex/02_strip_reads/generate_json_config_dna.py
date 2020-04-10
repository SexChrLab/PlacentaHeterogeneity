# In this script, I want to generate the json config file for stripping reads for dna

import json
import os.path
from os import path

data = {}

data['all_dna_samples'] = []

with open('c://Users/tuyen/Documents/postdoc_asu/projects/PlacentaSexDiff/B_gtex/01_download_data/wes_samples_females_final.csv', 'r') as f:
    for line in f:
        items = line.rstrip('\n').split(',')
        for i in items[1:]:
            data['all_dna_samples'].append(i)

data['GTEX_bam_dir'] = '/data/mwilsons/public_data/controlled_access/gtex/version8/WES_BAM'

with open(
    'c://Users/tuyen/Documents/postdoc_asu/projects/PlacentaSexDiff/B_gtex/02_strip_reads/strip_dna_read_config.json',
    'w') as outfile:
    json.dump(data, outfile)
