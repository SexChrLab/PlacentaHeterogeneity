import os

configfile: "asereadcounter_config.json"

chromosomes = ["8", "X"]

gatk_path = "/home/tphung3/softwares/gatk-4.1.0.0/gatk"

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
    input:
        expand("asereadcounter/HISAT/before_vqsr/chr{chr}/{combo}_chr{chr}.tsv", combo=combiList_placenta, chr=chromosomes)

rule gatk_asereadcounter_placenta:
    input:
        ref = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa",
        bam = "/data/storage/SAYRES/placenta_YPOPS/06_HISAT_RNA/Aligned_BAMS/stranded_RF/{rna}_HISAT_pair_trim_sort_mkdup_rdgrp_XX.bam",
        sites = "/scratch/tphung3/PlacentaSexDiff/A_placenta/05_genotype_rnaseq/genotyped_vcfs/chr{chr}.gatk.called.raw.biallelic.snp.{dna}.het.vcf"
    output:
        "asereadcounter/HISAT/before_vqsr/chr{chr}/{dna}_{rna}_chr{chr}.tsv"
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
