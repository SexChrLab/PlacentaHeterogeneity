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

### 04_process_rna

### 05_process_dna

### 06_run_asereadcounter

### 07_analyze_ase
