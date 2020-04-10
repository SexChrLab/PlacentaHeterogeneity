import os

configfile: "analyze_ase_config.json"

chromosomes = ["8", "X"]

import itertools

# placenta
combiC_placenta = []
for key in config["dna_rna_placenta"]:
    for item in config["dna_rna_placenta"][key]:
        combiC_placenta.append((key, item))

combiList_placenta=list()
for c in combiC_placenta:
    combiList_placenta.append(c[0]+"_"+c[1])

# decidua females
combiC_decidua_females = []
for key in config["dna_rna_decidua_females"]:
    for item in config["dna_rna_decidua_females"][key]:
        combiC_decidua_females.append((key, item))

combiList_decidua_females=list()
for c in combiC_decidua_females:
    combiList_decidua_females.append(c[0]+"_"+c[1])

# decidua males
combiC_decidua_males = []
for key in config["dna_rna_decidua_males"]:
    for item in config["dna_rna_decidua_males"][key]:
        combiC_decidua_males.append((key, item))

combiList_decidua_males=list()
for c in combiC_decidua_males:
    combiList_decidua_males.append(c[0]+"_"+c[1])

rule all:
    input: #calc_allele_balance
        expand("/scratch/tphung3/PlacentaSexDiff/A_placenta/03_analyze_ase/results/chr{chr}/{combo}_chr{chr}_allele_balance.tsv", combo=combiList_placenta, chr=chromosomes),
        expand("/scratch/tphung3/PlacentaSexDiff/A_placenta/03_analyze_ase/results/chr{chr}/{combo}_chr{chr}_allele_balance.tsv", combo=combiList_decidua_females, chr=chromosomes),
        expand("/scratch/tphung3/PlacentaSexDiff/A_placenta/03_analyze_ase/results/chr{chr}/{combo}_chr{chr}_allele_balance.tsv", combo=combiList_decidua_males, chr=chromosomes)

rule calc_allele_balance:
    input:
        "/scratch/tphung3/PlacentaSexDiff/A_placenta/02_run_asereadcounter/asereadcounter/HISAT/chr{chr}/{dna}_{rna}_chr{chr}.tsv"
    output:
        "/scratch/tphung3/PlacentaSexDiff/A_placenta/03_analyze_ase/results/chr{chr}/{dna}_{rna}_chr{chr}_allele_balance.tsv"
    params:
        script = config["calc_allele_balance_script"]
    shell:
        """
        python {params.script} --input {input} --output {output}
        """
