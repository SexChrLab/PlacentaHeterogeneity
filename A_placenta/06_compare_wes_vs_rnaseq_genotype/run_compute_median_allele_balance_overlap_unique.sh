#!/bin/bash

python /scratch/tphung3/PlacentaSexDiff/A_placenta/06_compare_wes_vs_rnaseq_genotype/scripts/compute_median_allele_balance_overlap_unique.py --WES_FILES_DIR /scratch/tphung3/PlacentaSexDiff/A_placenta/03_analyze_ase/results/chr8 --RNA_FILES_DIR /scratch/tphung3/PlacentaSexDiff/A_placenta/05_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chr8/ --threshold 10 --include_pars yes --outfile /scratch/tphung3/PlacentaSexDiff/A_placenta/06_compare_wes_vs_rnaseq_genotype/results/compare_wes_rna_genotyping_chr8.txt

awk 'NF' /scratch/tphung3/PlacentaSexDiff/A_placenta/06_compare_wes_vs_rnaseq_genotype/results/compare_wes_rna_genotyping_chr8.txt > /scratch/tphung3/PlacentaSexDiff/A_placenta/06_compare_wes_vs_rnaseq_genotype/results/compare_wes_rna_genotyping_chr8_fmt.txt

python /scratch/tphung3/PlacentaSexDiff/A_placenta/06_compare_wes_vs_rnaseq_genotype/scripts/compute_median_allele_balance_overlap_unique.py --WES_FILES_DIR /scratch/tphung3/PlacentaSexDiff/A_placenta/03_analyze_ase/results/chrX --RNA_FILES_DIR /scratch/tphung3/PlacentaSexDiff/A_placenta/05_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chrX/ --threshold 10 --include_pars no --outfile /scratch/tphung3/PlacentaSexDiff/A_placenta/06_compare_wes_vs_rnaseq_genotype/results/compare_wes_rna_genotyping_chrX.txt

awk 'NF' /scratch/tphung3/PlacentaSexDiff/A_placenta/06_compare_wes_vs_rnaseq_genotype/results/compare_wes_rna_genotyping_chrX.txt > /scratch/tphung3/PlacentaSexDiff/A_placenta/06_compare_wes_vs_rnaseq_genotype/results/compare_wes_rna_genotyping_chrX_fmt.txt
