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
    input:
        expand("asereadcounter/HISAT/chr{chr}/{combo}_chr{chr}.tsv", combo=combiList_placenta, chr=chromosomes),
        expand("asereadcounter/HISAT/chr{chr}/{combo}_chr{chr}.tsv", combo=combiList_decidua_females, chr=chromosomes),
        expand("asereadcounter/HISAT/chr{chr}/{combo}_chr{chr}.tsv", combo=combiList_decidua_males, chr=chromosomes)

rule gatk_asereadcounter_placenta:
    input:
        ref = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa",
        bam = "/data/storage/SAYRES/placenta_YPOPS/06_HISAT_RNA/Aligned_BAMS/stranded_RF/{rna}_HISAT_pair_trim_sort_mkdup_rdgrp_XX.bam",
        sites = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/vqsr/chr{chr}.gatk.called.vqsr.sv.biallelic.snp.{dna}.het.vcf"
    output:
        "asereadcounter/HISAT/chr{chr}/{dna}_{rna}_chr{chr}.tsv"
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

rule gatk_asereadcounter_decidua:
    input:
        ref = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa",
        bam = "/data/storage/SAYRES/placenta_YPOPS/06_HISAT_RNA/Aligned_BAMS/stranded_RF/{rna}_HISAT_pair_trim_sort_mkdup_rdgrp_DEC.bam",
        sites = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/vqsr/chr{chr}.gatk.called.vqsr.sv.biallelic.snp.{dna}.het.vcf"
    output:
        "asereadcounter/HISAT/chr{chr}/{dna}_{rna}_chr{chr}.tsv"
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
