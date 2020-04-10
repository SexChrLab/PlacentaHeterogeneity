# In this script, we want to calculate median allele balance across all variants for each individual within placenta and within decidua
# Allow for threshold on total count
# Input file is "/scratch/tphung3/PlacentaSexDiff/A_placenta/03_analyze_ase/results/chr{chr}/{combo}_chr{chr}_allele_balance.tsv"

import os
import pandas as pd
import argparse
import json

parser = argparse.ArgumentParser(description='Calculate median allele balance for within placenta and within decidua.')
parser.add_argument('--chromosome',required=True,help='Input the chromosome number')
parser.add_argument('--files_dir',required=True,help='Input the path to the directory where the tsv file of the allele balance is. For example: /scratch/tphung3/PlacentaSexDiff/A_placenta/03_analyze_ase/results/chrX')
parser.add_argument('--threshold',required=True,help='Input the threshold to subset the total count with')
parser.add_argument('--include_pars',required=True,help='Yes if PARs are included in the calculation. No otherwise.')
parser.add_argument('--out_placenta',required=True,help='Input the path to the output file for the placenta tissues')
parser.add_argument('--out_decidua_females',required=True,help='Input the path to the output file for the decidua tissues where the fetus is females')
parser.add_argument('--out_decidua_males',required=True,help='Input the path to the output file for the decidua tissues where the fetus is males')

args = parser.parse_args()

chromosome = args.chromosome
FILES_DIR = args.files_dir
threshold = float(args.threshold)
out_placenta = args.out_placenta
out_decidua_females = args.out_decidua_females
out_decidua_males = args.out_decidua_males

if not os.path.exists(out_placenta):
    header = ['tissue', 'chr', 'rna_id', 'num_var', 'median_allele_balance', 'num_var_subset', 'median_allele_balance_subset']
    outfile_placenta = open(out_placenta, 'w')
    print ('\t'.join(header), file=outfile_placenta)
else:
    outfile_placenta = open(out_placenta, 'a')

if not os.path.exists(out_decidua_females):
    header = ['tissue', 'chr', 'rna_id', 'num_var', 'median_allele_balance', 'num_var_subset', 'median_allele_balance_subset']
    out_decidua_females = open(out_decidua_females, 'w')
    print ('\t'.join(header), file=out_decidua_females)
else:
    out_decidua_females = open(out_decidua_females, 'a')

if not os.path.exists(out_decidua_males):
    header = ['tissue', 'chr', 'rna_id', 'num_var', 'median_allele_balance', 'num_var_subset', 'median_allele_balance_subset']
    out_decidua_males = open(out_decidua_males, 'w')
    print ('\t'.join(header), file=out_decidua_males)
else:
    out_decidua_males = open(out_decidua_males, 'a')

placenta_rna_ids = []
decidua_females_rna_ids = []
decidua_males_rna_ids = []
with open('/scratch/tphung3/PlacentaSexDiff/A_placenta/03_analyze_ase/scripts/analyze_ase_config.json') as json_file:
    data = json.load(json_file)
    for placenta_dna_id in data['dna_rna_placenta']:
        for placenta_rna_id in data['dna_rna_placenta'][placenta_dna_id]:
            placenta_rna_ids.append(placenta_rna_id)
    for decidua_dna_id in data['dna_rna_decidua_females']:
        for decidua_rna_id in data['dna_rna_decidua_females'][decidua_dna_id]:
            decidua_females_rna_ids.append(decidua_rna_id)
    for decidua_dna_id in data['dna_rna_decidua_males']:
        for decidua_rna_id in data['dna_rna_decidua_males'][decidua_dna_id]:
            decidua_males_rna_ids.append(decidua_rna_id)

for file in [f for f in os.listdir(FILES_DIR) if f.endswith('_allele_balance.tsv')]:
    items = file.split('_')
    if items[1] in placenta_rna_ids:
        data = pd.read_csv(os.path.join(FILES_DIR, file), sep='\t') #header of input file is: header = ['chr', 'position', 'ref_allele', 'alt_allele', 'ref_count', 'alt_count', 'total_count', 'allele_balance']
        if args.include_pars == 'yes':
            data_nrow = len(data.index)
            median_allele_balance = data['allele_balance'].median()
            # Subset based on threshold of total count
            data_subset = data[data['total_count'] > threshold]
            data_subset_nrow = len(data_subset.index)
            median_allele_balance_subset = data_subset['allele_balance'].median()

            output = ['placenta', chromosome, items[1], str(data_nrow), str(median_allele_balance), str(data_subset_nrow), str(median_allele_balance_subset)]
            print ('\t'.join(output), file=outfile_placenta)
        else:  # remove pars
            data_rmpars_1 = data[(data['position'] < 10001) | (data['position'] > 2781479)]
            data_rmpars_both = data_rmpars_1[
                (data_rmpars_1['position'] < 155701383) | (data_rmpars_1['position'] > 156030895)]
            data_rmpars_both_nrow = len(data_rmpars_both.index)
            # print (str(data_nrow))
            median_allele_balance = data_rmpars_both['allele_balance'].median()
            # Subset based on threshold of total count
            data_rmpars_both_subset = data_rmpars_both[data_rmpars_both['total_count'] > threshold]
            data_rmpars_both_subset_nrow = len(data_rmpars_both_subset.index)
            median_allele_balance_subset = data_rmpars_both_subset['allele_balance'].median()

            output = ['placenta', chromosome, items[1], str(data_rmpars_both_nrow), str(median_allele_balance),
                      str(data_rmpars_both_subset_nrow), str(median_allele_balance_subset)]
            print('\t'.join(output), file=outfile_placenta)
    elif items[1] in decidua_females_rna_ids:
        if args.include_pars == 'yes':
            data_nrow = len(data.index)
            median_allele_balance = data['allele_balance'].median()
            # Subset based on threshold of total count
            data_subset = data[data['total_count'] > threshold]
            data_subset_nrow = len(data_subset.index)
            median_allele_balance_subset = data_subset['allele_balance'].median()

            output = ['decidua', chromosome, items[1], str(data_nrow), str(median_allele_balance), str(data_subset_nrow), str(median_allele_balance_subset)]
            print ('\t'.join(output), file=out_decidua_females)
        else:  # remove pars
            data_rmpars_1 = data[(data['position'] < 10001) | (data['position'] > 2781479)]
            data_rmpars_both = data_rmpars_1[
                (data_rmpars_1['position'] < 155701383) | (data_rmpars_1['position'] > 156030895)]
            data_rmpars_both_nrow = len(data_rmpars_both.index)
            # print (str(data_nrow))
            median_allele_balance = data_rmpars_both['allele_balance'].median()
            # Subset based on threshold of total count
            data_rmpars_both_subset = data_rmpars_both[data_rmpars_both['total_count'] > threshold]
            data_rmpars_both_subset_nrow = len(data_rmpars_both_subset.index)
            median_allele_balance_subset = data_rmpars_both_subset['allele_balance'].median()

            output = ['decidua', chromosome, items[1], str(data_rmpars_both_nrow), str(median_allele_balance),
                      str(data_rmpars_both_subset_nrow), str(median_allele_balance_subset)]
            print('\t'.join(output), file=out_decidua_females)
    else:
        if args.include_pars == 'yes':
            data_nrow = len(data.index)
            median_allele_balance = data['allele_balance'].median()
            # Subset based on threshold of total count
            data_subset = data[data['total_count'] > threshold]
            data_subset_nrow = len(data_subset.index)
            median_allele_balance_subset = data_subset['allele_balance'].median()

            output = ['decidua', chromosome, items[1], str(data_nrow), str(median_allele_balance), str(data_subset_nrow), str(median_allele_balance_subset)]
            print ('\t'.join(output), file=out_decidua_males)
        else:  # remove pars
            data_rmpars_1 = data[(data['position'] < 10001) | (data['position'] > 2781479)]
            data_rmpars_both = data_rmpars_1[
                (data_rmpars_1['position'] < 155701383) | (data_rmpars_1['position'] > 156030895)]
            data_rmpars_both_nrow = len(data_rmpars_both.index)
            # print (str(data_nrow))
            median_allele_balance = data_rmpars_both['allele_balance'].median()
            # Subset based on threshold of total count
            data_rmpars_both_subset = data_rmpars_both[data_rmpars_both['total_count'] > threshold]
            data_rmpars_both_subset_nrow = len(data_rmpars_both_subset.index)
            median_allele_balance_subset = data_rmpars_both_subset['allele_balance'].median()

            output = ['decidua', chromosome, items[1], str(data_rmpars_both_nrow), str(median_allele_balance),
                      str(data_rmpars_both_subset_nrow), str(median_allele_balance_subset)]
            print('\t'.join(output), file=out_decidua_males)

