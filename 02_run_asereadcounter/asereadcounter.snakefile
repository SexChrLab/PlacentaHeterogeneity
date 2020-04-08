import os

configfile: "asereadcounter_config.json"

chromosomes = ["8", "X"]

gatk_path = "/home/tphung3/softwares/gatk-4.1.0.0/gatk"

combiC_placenta = []
for key in config["dna_rna_placenta"]:
    for item in config["dna_rna_placenta"][key]:
        combiC_placenta.append((key, item))

import itertools
combiList_placenta=list()
for c in combiC_placenta:
    combiList_placenta.append(c[0]+"_"+c[1])

combiC_decidua = []
for key in config["dna_rna_decidua"]:
    for item in config["dna_rna_decidua"][key]:
        combiC_decidua.append((key, item))

import itertools
combiList_decidua=list()
for c in combiC_decidua:
    combiList_decidua.append(c[0]+"_"+c[1])

rule all:
    input: #Run asereadcounter for autosomes and chrX for decidua samples
        expand("asereadcounter/HISAT/chr{chr}/{combo}_chr{chr}.tsv", combo=combiList_decidua, chr=chromosomes)
    # input: #Run asereadcounter for autosomes and chrX for placenta samples
    #     expand("asereadcounter/HISAT/chr{chr}/{combo}_chr{chr}.tsv", combo=combiList_placenta, chr=chromosomes)

# rule gatk_asereadcounter_placenta:
#     input:
#         ref = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa",
#         bam = "/data/storage/SAYRES/placenta_YPOPS/06_HISAT_RNA/Aligned_BAMS/stranded_RF/{rna}_HISAT_pair_trim_sort_mkdup_rdgrp_XX.bam",
#         sites = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/vqsr/chr{chr}.gatk.called.vqsr.sv.biallelic.snp.{dna}.het.vcf"
#     output:
#         "asereadcounter/HISAT/chr{chr}/{dna}_{rna}_chr{chr}.tsv"
#     params:
#         gatk = gatk_path
#     shell:
#         """{params.gatk} ASEReadCounter """
#         """-R {input.ref} """
#         """--output {output} """
#         """--input {input.bam} """
#         """--variant {input.sites} """
#         """--min-depth-of-non-filtered-base 1 """
#         """--min-mapping-quality 10 """
#         """--min-base-quality 10 """

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
