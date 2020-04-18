import json
from collections import defaultdict

new_data = {}
new_data['dna_rna'] = defaultdict(list)

with open('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/genotype_rnaseq_config.json') as json_file:
    data = json.load(json_file)
    for i in data['all_rna_samples']:
        rna = i
        dna = i.split('-SM-')[0]
        new_data['dna_rna'][dna].append(rna)

with open('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/asereadcounter_config.json', 'w') as outfile:
    json.dump(new_data, outfile)
