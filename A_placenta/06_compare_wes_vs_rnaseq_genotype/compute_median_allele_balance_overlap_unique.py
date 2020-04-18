# In this script, we want to calculate the median allele balance for variants that are shared between the 2 sets (genotyping using wes vs genotyping using rnaseq) and for variants that are unique to each set.

import os
import sys
import pandas as pd
import argparse



parser = argparse.ArgumentParser(description='Compute median allele balance to compare between genotyping from wes versus rnaseq.')
parser.add_argument('--WES_FILES_DIR',required=True,help='Input the path to the asereadcounter file after calculating allele balance for genotyping from wes. For example: ~/scratch/PlacentaSexDiff/A_placenta/03_analyze_ase/results/chrX/')
parser.add_argument('--RNA_FILES_DIR',required=True,help='Input the path to the asereadcounter file after calculating allele balance for genotyping from rnaseq. For example: ~/scratch/PlacentaSexDiff/A_placenta/05_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chrX/')
parser.add_argument('--threshold',required=True,help='Input the threshold to subset the total count with')
parser.add_argument('--include_pars',required=True,help='Yes if PARs are included in the calculation. No otherwise.')
parser.add_argument('--outfile',required=True,help='Input the path to the output file')

args = parser.parse_args()

WES_FILES_DIR = args.WES_FILES_DIR
RNA_FILES_DIR = args.RNA_FILES_DIR
threshold = float(args.threshold)

# Set up the output file
if not os.path.exists(args.outfile):
    header = ['wes_fn', 'rna_fn', 'wes_numvar', 'rna_numvar', 'wes_subset_numvar', 'rna_subset_numvar', 'wes_subset_median', 'rna_subset_median', 'shared_numvar', 'wes_subset_shared_median', 'rna_subset_shared_median', 'wes_subset_unique_numvar', 'rna_subset_unique_numvar', 'wes_subset_unique_median', 'rna_subset_unique_median']
    outfile = open(args.outfile, 'w')
    print ('\t'.join(header), file=outfile)
else:
    outfile = open(args.outfile, 'a')


for wes_file in [f for f in os.listdir(WES_FILES_DIR) if f.endswith('_allele_balance.tsv')]: #OBG0111_OBG0111-2-011_chr8_allele_balance.tsv from wes vs OBG0111-2-011_RNA_OBG0111-2-011_chr8_allele_balance.tsv from rnaseq
    out = []
    rna_name = wes_file.split('_')[1]

    for rna_file in [f for f in os.listdir(RNA_FILES_DIR) if f.endswith('_allele_balance.tsv')]:
        if rna_file.startswith(rna_name):
            out.append(wes_file)
            out.append(rna_file)
            # Open file in a panda dataframe
            wes_data = pd.read_csv(os.path.join(WES_FILES_DIR, wes_file),
                                   sep='\t')  # header of input file is: header = ['chr', 'position', 'ref_allele', 'alt_allele', 'ref_count', 'alt_count', 'total_count', 'allele_balance']
            rna_data = pd.read_csv(os.path.join(RNA_FILES_DIR, rna_file), sep='\t')

            if args.include_pars == 'yes':
                out.append(str(wes_data.shape[0]))
                out.append(str(rna_data.shape[0]))

                # Filter based on threshold
                wes_data_subset = wes_data[wes_data['total_count'] > 10]
                out.append(str(wes_data_subset.shape[0]))

                rna_data_subset = rna_data[rna_data['total_count'] > 10]
                out.append(str(rna_data_subset.shape[0]))

                # Compute median allele balance for the subset data
                out.append(str(wes_data_subset['allele_balance'].median()))
                out.append(str(rna_data_subset['allele_balance'].median()))

                # Find shared variants
                common_variants = pd.merge(wes_data_subset, rna_data_subset, on='position')
                out.append(str(common_variants.shape[0]))

                out.append(str(common_variants['allele_balance_x'].median()))
                out.append(str(common_variants['allele_balance_y'].median()))

                # Find unique variants to each set
                wes_data_unique = wes_data_subset[~wes_data_subset['position'].isin(common_variants['position'])].dropna()
                rna_data_unique = rna_data_subset[~rna_data_subset['position'].isin(common_variants['position'])].dropna()

                out.append(str(wes_data_unique.shape[0]))
                out.append(str(rna_data_unique.shape[0]))

                out.append(str(wes_data_unique['allele_balance'].median()))
                out.append(str(rna_data_unique['allele_balance'].median()))
            else:
                wes_data_rmpars_1 = wes_data[(wes_data['position'] < 10001) | (wes_data['position'] > 2781479)]
                wes_data_rmpars_both = wes_data_rmpars_1[
                    (wes_data_rmpars_1['position'] < 155701383) | (wes_data_rmpars_1['position'] > 156030895)]

                rna_data_rmpars_1 = rna_data[(rna_data['position'] < 10001) | (rna_data['position'] > 2781479)]
                rna_data_rmpars_both = rna_data_rmpars_1[
                    (rna_data_rmpars_1['position'] < 155701383) | (rna_data_rmpars_1['position'] > 156030895)]

                out.append(str(wes_data_rmpars_both.shape[0]))
                out.append(str(rna_data_rmpars_both.shape[0]))

                # Filter based on threshold
                wes_data_subset = wes_data_rmpars_both[wes_data_rmpars_both['total_count'] > 10]
                out.append(str(wes_data_subset.shape[0]))

                rna_data_subset = rna_data_rmpars_both[rna_data_rmpars_both['total_count'] > 10]
                out.append(str(rna_data_subset.shape[0]))

                # Compute median allele balance for the subset data
                out.append(str(wes_data_subset['allele_balance'].median()))
                out.append(str(rna_data_subset['allele_balance'].median()))

                # Find shared variants
                common_variants = pd.merge(wes_data_subset, rna_data_subset, on='position')
                out.append(str(common_variants.shape[0]))

                out.append(str(common_variants['allele_balance_x'].median()))
                out.append(str(common_variants['allele_balance_y'].median()))

                # Find unique variants to each set
                wes_data_unique = wes_data_subset[~wes_data_subset['position'].isin(common_variants['position'])].dropna()
                rna_data_unique = rna_data_subset[~rna_data_subset['position'].isin(common_variants['position'])].dropna()

                out.append(str(wes_data_unique.shape[0]))
                out.append(str(rna_data_unique.shape[0]))

                out.append(str(wes_data_unique['allele_balance'].median()))
                out.append(str(rna_data_unique['allele_balance'].median()))


    print ('\t'.join(out), file=outfile)


