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
