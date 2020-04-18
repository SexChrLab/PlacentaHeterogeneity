# PlacentaHeterogeneity

## A_placenta
### 01_data_processing
**Whole exome**
- Config file: `process_dna_females_config.json`
  - The config file is generated with `python generate_json_config_dna_females.py`. The csv file `female_sample_ids.csv` is needed for the python script `generate_json_config_dna_females.py`: the first column is the sample id, the second column is the name of the fastq file for read 1, and the third column is the name of the fastq file for read 2. 
- Snakefile: `process_dna_females.snakefile`
  
### 02_run_asereadcounter
- Snakefile: `asereadcounter.snakefile`
- Config file: `asereadcounter_config.json`

### 03_analyze_ase
1. Calculate allele balance (unphased)
- `calc_allele_balance.py`
- Snakefile `analyze_ase.snakefile` (config file `analyze_ase_config.json`)

2. Calculate median allele balance for each individual
- Remove any variant where the total count is less than 10
- `calc_median_allele_balance_placenta_decidua.py`
- Command line: `run_calc_median_allele_balance_placenta_decidua.sh`

### 04_analyze_xist
- From the results of running GATK ASEReadcounter, we want to select variants that are in the XIST region (chrX:73820892-73851867)
- `subset_for_positions_in_xist.py`: This script selects all of the outputs from GATK ASEReadCounter and returns any variant within the XIST region. 

### 05_genotype_rnaseq
- Genotype calling using the RNAseq data
- Placenta samples only for now 
1. Snakefile `genotype_rnaseq.snakefile` (config file `genotype_rnaseq_config.json`): use GATK joint genotyping
- Note that the filtering with VQSR does not work. For now, we are working with the set of raw variants
2. Snakefile `asereadcounter.snakefile` (config file `asereadcounter_config.json`): run asereadcounter
3. Snakefile `analyze_ase.snakefile` (config file `asereadcounter_config.json`): calculate allele balance
4. Python script `calc_median_allele_balance_placenta.py` (command line `run_calc_median_allele_balance_placenta_decidua.sh`): calculate median allele balance for each individual

### 06_compare_wes_vs_rnaseq_genotype
- Compare heterozygous and expressed variants identified from using the whole exome data for genotyping versus from using the rnaseq data for genotyping
- Use the Python script `compute_median_allele_balance_overlap_unique.py` (see bash script `run_compute_median_allele_balance_overlap_unique.sh`)
- Plotting using R script `compare_wes_rnaseq.R`

## B_gtex
### 01_download_data
1. Find females in GTEx: `find_females_in_GTEx.py`
2. Parse GTEx information to subset to different tissues: `parse_GTEx_info.py`
3. RNA sample ids: `rna_samples_females_final.csv`
4. DNA sample ids: `wes_samples_females_final`

### 02_strip_reads
1. Strip RNAseq reads
  1. Use the Python script `generate_json_config_rna.py`
  2. Use snakefile `strip_rna_read.snakefile`
2. Strip WES reads
  1. Use the Python script `generate_json_config_dna.py`
  2. Use the snakefile `strip_dna_read.snakefile`
3. QC: unequal number of reads in RNAseq data
- In some RNAseq samples, the number of R1 reads and R2 reads are unequal. Therefore, I wrote a Python script `count_number_of_reads.py` to count the number of reads
- The script is implemented in the snakefile `count_number_of_reads.snakefile`
- Find the samples with unequal number of reads in order to remove them:
  - `return_samples_unequal_number_of_reads.py`

### 03_extract_RG
1. `cp PlacentaSexDiff/B_gtex/02_strip_reads/strip_rna_read_config.json PlacentaSexDiff/B_gtex/04_process_rna/process_rna_config.json`
2. `cp PlacentaSexDiff/B_gtex/02_strip_reads/strip_dna_read_config.json PlacentaSexDiff/B_gtex/05_process_dna/process_dna_config.json`
3. Extract read groups:
  - `snakemake --snakefile extract_RG_rna.snakefile`
  - `snakemake --snakefile extract_RG_dna.snakefile`
4. Add to json files for the 04_process_rna step and 05_process_dna step
  - `python make_json_for_read_groups_rna.py`
  - `python make_json_for_read_groups_dna.py`

### 04_process_rna
- Snakefile: `process_rna.snakefile`
- Config file: `process_rna_config.json`

### 05_process_dna
- Snakefile: `process_dna.snakefile`
- Config file: `process_dna_config.json`
- After genotyping, I want to use bcftools in order to subset each VCF file for each individual for running asereadcounter.
  - However, the problem is that the ID in the list `all_dna_samples` is like this: `GTEX-111CU-0003-SM-58Q95` while the header for the VCF file is like this `GTEX-111CU-0003`.
  - Therefore, I created another list called `all_dna_samples_without_SM` to contain modified ID `GTEX-111CU-0003`. Use the script `modify_json.py`

### 06_run_asereadcounter
- Generate the config file with script `make_json.py`
- Snakefile: `asereadcounter.snakefile`
- Config file: `asereadcounter_config.json`

### 07_analyze_ase
- Config file: `analyze_ase_config.json`
- Snakefile: `analyze_ase.snakefile`
  - Calculate allele balance (unphased): `calc_allele_balance.py`
    - rule `calc_allele_balance`
  - Calculate median allele balance for each individual per tissue
    - Remove any variant where the total count is less than 10
    - `calc_median_allele_balance_placenta_decidua.py`
    - rule `calc_median_allele_balance_per_tissue_chr8` and `calc_median_allele_balance_per_tissue_chrX_nonpars`
