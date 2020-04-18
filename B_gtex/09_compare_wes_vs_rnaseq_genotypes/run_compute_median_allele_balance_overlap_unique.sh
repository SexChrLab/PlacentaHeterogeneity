#!/bin/bash

for tissue in AdrenalGland ArteryAorta BrainAmygdala BrainCortex BreastMammaryTissue ColonSigmoid HeartLeftVentricle KidneyCortex Liver Lung MuscleSkeletal NerveTibial Ovary Pancreas Pituitary Spleen Stomach Thyroid Uterus Vagina
do

python /scratch/tphung3/PlacentaSexDiff/B_gtex/09_compare_wes_vs_rnaseq_genotypes/scripts/compute_median_allele_balance_overlap_unique.py --WES_FILES_DIR /scratch/tphung3/PlacentaSexDiff/B_gtex/07_analyze_ase/results/chr8/ --RNA_FILES_DIR /scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chr8/ --threshold 10 --include_pars yes --tissue ${tissue} --outfile /scratch/tphung3/PlacentaSexDiff/B_gtex/09_compare_wes_vs_rnaseq_genotypes/results/chr8/compare_wes_rna_genotyping_chr8_${tissue}.txt

python /scratch/tphung3/PlacentaSexDiff/B_gtex/09_compare_wes_vs_rnaseq_genotypes/scripts/compute_median_allele_balance_overlap_unique.py --WES_FILES_DIR /scratch/tphung3/PlacentaSexDiff/B_gtex/07_analyze_ase/results/chrX/ --RNA_FILES_DIR /scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chrX/ --threshold 10 --include_pars no --tissue ${tissue} --outfile /scratch/tphung3/PlacentaSexDiff/B_gtex/09_compare_wes_vs_rnaseq_genotypes/results/chrX/compare_wes_rna_genotyping_chrX_${tissue}.txt

done
