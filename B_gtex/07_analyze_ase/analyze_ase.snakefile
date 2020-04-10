import os

configfile: "analyze_ase_config.json"

chromosomes = ["8", "X"]

combiC = []
for key in config["dna_rna"]:
    for item in config["dna_rna"][key]:
        combiC.append((key, item))

import itertools
combiList=list()
for c in combiC:
    combiList.append(c[0]+"_"+c[1])

rule all:
    input: #calc_median_allele_balance_per_tissue
        expand("/scratch/tphung3/PlacentaSexDiff/B_gtex/07_analyze_ase/results/chr8/{tissue}_chr8_median_allele_balance.txt", tissue=config["tissues"]), #chr8
        expand("/scratch/tphung3/PlacentaSexDiff/B_gtex/07_analyze_ase/results/chrX/{tissue}_chrX_nonpars_median_allele_balance.txt", tissue=config["tissues"])
    input: #calc_allele_balance
        expand("/scratch/tphung3/PlacentaSexDiff/B_gtex/07_analyze_ase/results/chr{chr}/{combo}_chr{chr}_allele_balance.tsv", combo=combiList, chr=chromosomes)

rule calc_allele_balance:
    input:
        "/scratch/tphung3/PlacentaSexDiff/B_gtex/06_run_asereadcounter/asereadcounter/HISAT/chr{chr}/{dna}_{rna}_chr{chr}.tsv"
    output:
        "/scratch/tphung3/PlacentaSexDiff/B_gtex/07_analyze_ase/results/chr{chr}/{dna}_{rna}_chr{chr}_allele_balance.tsv"
    params:
        script = config["calc_allele_balance_script"]
    shell:
        """
        python {params.script} --input {input} --output {output}
        """

rule calc_median_allele_balance_per_tissue_chr8:
    output:
        "/scratch/tphung3/PlacentaSexDiff/B_gtex/07_analyze_ase/results/chr8/{tissue}_chr8_median_allele_balance.txt"
    params:
        script = config["calc_median_allele_balance_per_tissue_script"],
        tissue = "{tissue}",
        chr = "chr8",
        dir = "/scratch/tphung3/PlacentaSexDiff/B_gtex/07_analyze_ase/results/chr8",
        threshold = 10
    shell:
        """
        python {params.script} --tissue {params.tissue} --chromosome {params.chr} --files_dir {params.dir} --threshold {params.threshold} --include_pars yes --out {output}
        """

rule calc_median_allele_balance_per_tissue_chrX_nonpars:
    output:
        "/scratch/tphung3/PlacentaSexDiff/B_gtex/07_analyze_ase/results/chrX/{tissue}_chrX_nonpars_median_allele_balance.txt"
    params:
        script = config["calc_median_allele_balance_per_tissue_script"],
        tissue = "{tissue}",
        chr = "chrX",
        dir = "/scratch/tphung3/PlacentaSexDiff/B_gtex/07_analyze_ase/results/chrX",
        threshold = 10
    shell:
        """
        python {params.script} --tissue {params.tissue} --chromosome {params.chr} --files_dir {params.dir} --threshold {params.threshold} --include_pars no --out {output}
        """
