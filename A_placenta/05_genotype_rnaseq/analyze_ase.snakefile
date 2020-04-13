import os

configfile: "asereadcounter_config.json"

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

rule all:
    input: #calc_allele_balance
        expand("asereadcounter/HISAT/before_vqsr/chr{chr}/{combo}_chr{chr}_allele_balance.tsv", combo=combiList_placenta, chr=chromosomes)

rule calc_allele_balance:
    input:
        "asereadcounter/HISAT/before_vqsr/chr{chr}/{combo}_chr{chr}.tsv"
    output:
        "asereadcounter/HISAT/before_vqsr/chr{chr}/{combo}_chr{chr}_allele_balance.tsv"
    shell:
        """
        python /scratch/tphung3/PlacentaSexDiff/A_placenta/03_analyze_ase/scripts/calc_allele_balance.py --input {input} --output {output}
        """
