#!/bin/bash

python calc_median_allele_balance_placenta_decidua.py --chromosome chr8 --files_dir ../results/chr8 --threshold 10 --out_placenta ../results/chr8/placenta_chr8_median_allele_balance.txt --out_decidua_females ../results/chr8/decidua_females_chr8_median_allele_balance.txt --out_decidua_males ../results/chr8/decidua_males_chr8_median_allele_balance.txt --include_pars yes

python calc_median_allele_balance_placenta_decidua.py --chromosome chrX --files_dir ../results/chrX --threshold 10 --out_placenta ../results/chrX/placenta_chrX_nonpars_median_allele_balance.txt --out_decidua_females ../results/chrX/decidua_females_chrX_nonpars_median_allele_balance.txt --out_decidua_males ../results/chrX/decidua_males_chrX_nonpars_median_allele_balance.txt --include_pars no
