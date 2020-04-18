import os

configfile: "asereadcounter_config.json"

chromosomes = ["8", "X"]

gatk_path = "/home/tphung3/softwares/gatk-4.1.0.0/gatk"

import itertools

combiC = []
for key in config["dna_rna"]:
    for item in config["dna_rna"][key]:
        combiC.append((key, item))

combiList=list()
for c in combiC:
    combiList.append(c[0]+"_"+c[1])

rule all:
    input: #after vqsr
        expand("asereadcounter/HISAT/after_vqsr/chr{chr}/{combo}_chr{chr}.tsv", combo=combiList, chr=chromosomes)

rule gatk_asereadcounter_placenta_after_vqsr:
    input:
        ref = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa",
        bam = "/scratch/tphung3/PlacentaSexDiff/B_gtex/04_process_rna/processed_bams/rna/{rna}.GRCh38.p12.genome.XXonly.sorted.merged.bam",
        sites = "/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/vqsr/chr{chr}.gatk.called.vqsr.sv.biallelic.snp.{dna}.het.vcf"
    output:
        "asereadcounter/HISAT/after_vqsr/chr{chr}/{dna}_{rna}_chr{chr}.tsv"
    params:
        gatk = gatk_path
    shell:
        """{params.gatk} ASEReadCounter """
        """-R {input.ref} """
        """--output {output} """
        """--input {input.bam} """
        """--variant {input.sites} """
        """--min-depth-of-non-filtered-base 1 """
        """--min-mapping-quality 10 """
        """--min-base-quality 10 """
