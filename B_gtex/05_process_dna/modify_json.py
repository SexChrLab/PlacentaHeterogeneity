# In this script, I want to modify the json file. Please see the explanation in the file notes.md for why I am doing this.
import json
with open('process_dna_config.json') as json_file:
    data = json.load(json_file)
data['all_dna_samples_without_SM'] = []
for i in data['all_dna_samples']:
    items = i.split('-')
    data['all_dna_samples_without_SM'].append(items[0] + '-' + items[1] + '-' + items[2])

with open('process_dna_config.json', 'w') as json_file:
    json.dump(data, json_file)

