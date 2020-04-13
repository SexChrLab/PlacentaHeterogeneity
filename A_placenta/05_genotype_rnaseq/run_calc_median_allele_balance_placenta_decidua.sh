#!/bin/bash

python calc_median_allele_balance_placenta.py --chromosome chr8 --files_dir /scratch/tphung3/PlacentaSexDiff/A_placenta/05_genotype_rnaseq/asereadcounter/HISAT/before_vqsr/chr8 --threshold 10 --out_placenta /scratch/tphung3/PlacentaSexDiff/A_placenta/05_genotype_rnaseq/asereadcounter/HISAT/before_vqsr/chr8/placenta_chr8_median_allele_balance.txt --include_pars yes

python calc_median_allele_balance_placenta.py --chromosome chrX --files_dir /scratch/tphung3/PlacentaSexDiff/A_placenta/05_genotype_rnaseq/asereadcounter/HISAT/before_vqsr/chrX --threshold 10 --out_placenta /scratch/tphung3/PlacentaSexDiff/A_placenta/05_genotype_rnaseq/asereadcounter/HISAT/before_vqsr/chrX/placenta_chrX_median_allele_balance.txt --include_pars no
