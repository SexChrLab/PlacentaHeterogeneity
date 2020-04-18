# In this script, I am generating the config file for running the snakefile genotype_rnaseq.snakefile

import json

new_data = {}
new_data['all_rna_samples'] = []
new_data['all_rna_samples_short'] = []

with open('/scratch/tphung3/PlacentaSexDiff/B_gtex/04_process_rna/process_rna_config.json') as json_file:
    data = json.load(json_file)
    for i in data['all_rna_samples']:
        new_data['all_rna_samples'].append(i)
        items = i.split('-SM-')
        new_data['all_rna_samples_short'].append(items[0])


with open('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/genotype_rnaseq_config.json', 'w') as outfile:
    json.dump(new_data, outfile)