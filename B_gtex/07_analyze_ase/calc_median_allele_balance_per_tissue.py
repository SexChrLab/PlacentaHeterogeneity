# In this script, we want to calculate median allele balance across all variants for each individual per tissue
# Allow for threshold on total count
# Input file is "/scratch/tphung3/PlacentaSexDiff/B_gtex/07_analyze_ase/results/chr{chr}/{combo}_chr{chr}_allele_balance.tsv"

import sys
import os
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Calculate median allele balance per tissue.')
parser.add_argument('--tissue',required=True,help='Input the name of the tissue')
parser.add_argument('--chromosome',required=True,help='Input the chromosome number')
parser.add_argument('--files_dir',required=True,help='Input the path to the directory where the tsv file of the allele balance is. For example: /scratch/tphung3/PlacentaSexDiff/B_gtex/07_analyze_ase/results/chrX')
parser.add_argument('--threshold',required=True,help='Input the threshold to subset the total count with')
parser.add_argument('--include_pars',required=True,help='Yes if PARs are included in the calculation. No otherwise.')
parser.add_argument('--out',required=True,help='Input the path to the output file')

args = parser.parse_args()

tissue = args.tissue
chromosome = args.chromosome
FILES_DIR = args.files_dir
threshold = float(args.threshold)
out = args.out

if not os.path.exists(out):
    header = ['tissue', 'chr', 'rna_id', 'num_var', 'median_allele_balance', 'num_var_subset', 'median_allele_balance_subset']
    outfile = open(out, 'w')
    print ('\t'.join(header), file=outfile)
else:
    outfile = open(out, 'a')

with open('/scratch/tphung3/PlacentaSexDiff/B_gtex/01_download_data/rna_samples_females_final_forcode.csv', 'r') as f:
    for line in f:
        if line.startswith(tissue):
            rna_ids = line.rstrip('\n').split(',')[1:]

# for file in [os.path.join(FILES_DIR, f) for f in os.listdir(FILES_DIR) if f.endswith('.tsv')]:
#     print (file)

for file in [f for f in os.listdir(FILES_DIR) if f.endswith('_allele_balance.tsv')]:
    items = file.split('_')
    if items[1] in rna_ids:
        data = pd.read_csv(os.path.join(FILES_DIR, file), sep='\t') #header of input file is: header = ['chr', 'position', 'ref_allele', 'alt_allele', 'ref_count', 'alt_count', 'total_count', 'allele_balance']
        if args.include_pars == 'yes':
            data_nrow = len(data.index)
            # print (str(data_nrow))
            median_allele_balance = data['allele_balance'].median()
            # Subset based on threshold of total count
            data_subset = data[data['total_count'] > threshold]
            data_subset_nrow = len(data_subset.index)
            median_allele_balance_subset = data_subset['allele_balance'].median()

            output = [tissue, chromosome, items[1], str(data_nrow), str(median_allele_balance), str(data_subset_nrow), str(median_allele_balance_subset)]
            print ('\t'.join(output), file=outfile)
        else: #remove pars
            data_rmpars_1 = data[(data['position'] < 10001) | (data['position'] > 2781479)]
            data_rmpars_both = data_rmpars_1[(data_rmpars_1['position'] < 155701383) | (data_rmpars_1['position'] > 156030895)]
            data_rmpars_both_nrow = len(data_rmpars_both.index)
            # print (str(data_nrow))
            median_allele_balance = data_rmpars_both['allele_balance'].median()
            # Subset based on threshold of total count
            data_rmpars_both_subset = data_rmpars_both[data_rmpars_both['total_count'] > threshold]
            data_rmpars_both_subset_nrow = len(data_rmpars_both_subset.index)
            median_allele_balance_subset = data_rmpars_both_subset['allele_balance'].median()

            output = [tissue, chromosome, items[1], str(data_rmpars_both_nrow), str(median_allele_balance), str(data_rmpars_both_subset_nrow), str(median_allele_balance_subset)]
            print ('\t'.join(output), file=outfile)