import os

configfile: "genotype_rnaseq_config.json"

chr_to_genotype = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"]
chr_current = ["8", "X"]

gatk_path = "/home/tphung3/softwares/gatk-4.1.0.0/gatk"
bcftools_path = "bcftools"


rule all:
    input: #further processing vcfs before vqsr
        expand("genotyped_vcfs/chr{chr}.gatk.called.raw.biallelic.snp.{rna}.het.vcf", chr=chr_current, rna=config["all_rna_samples_forbcftools"])
    input: #further processing vcfs before vqsr
        expand("genotyped_vcfs/chr{chr}.gatk.called.raw.biallelic.snp.{rna}.vcf", chr=chr_current, rna=config["all_rna_samples_forbcftools"])
    input: #GenotypeGVCFs
        expand("genotyped_vcfs/GRCh38.p12.genome.XXonly.chr{chr}.gatk.called.raw.vcf.gz", chr=chr_to_genotype)
    input: #combine gVCFs
        "combined_gvcfs/GRCh38.p12.genome.XXonly.chr1.gatk.combinegvcf.g.vcf.gz",
        "combined_gvcfs/GRCh38.p12.genome.XXonly.chr2.gatk.combinegvcf.g.vcf.gz",
        "combined_gvcfs/GRCh38.p12.genome.XXonly.chr3.gatk.combinegvcf.g.vcf.gz",
        "combined_gvcfs/GRCh38.p12.genome.XXonly.chr4.gatk.combinegvcf.g.vcf.gz",
        "combined_gvcfs/GRCh38.p12.genome.XXonly.chr5.gatk.combinegvcf.g.vcf.gz",
        "combined_gvcfs/GRCh38.p12.genome.XXonly.chr6.gatk.combinegvcf.g.vcf.gz",
        "combined_gvcfs/GRCh38.p12.genome.XXonly.chr7.gatk.combinegvcf.g.vcf.gz",
        "combined_gvcfs/GRCh38.p12.genome.XXonly.chr8.gatk.combinegvcf.g.vcf.gz",
        "combined_gvcfs/GRCh38.p12.genome.XXonly.chr9.gatk.combinegvcf.g.vcf.gz",
        "combined_gvcfs/GRCh38.p12.genome.XXonly.chr10.gatk.combinegvcf.g.vcf.gz",
        "combined_gvcfs/GRCh38.p12.genome.XXonly.chr11.gatk.combinegvcf.g.vcf.gz",
        "combined_gvcfs/GRCh38.p12.genome.XXonly.chr12.gatk.combinegvcf.g.vcf.gz",
        "combined_gvcfs/GRCh38.p12.genome.XXonly.chr13.gatk.combinegvcf.g.vcf.gz",
        "combined_gvcfs/GRCh38.p12.genome.XXonly.chr14.gatk.combinegvcf.g.vcf.gz",
        "combined_gvcfs/GRCh38.p12.genome.XXonly.chr15.gatk.combinegvcf.g.vcf.gz",
        "combined_gvcfs/GRCh38.p12.genome.XXonly.chr16.gatk.combinegvcf.g.vcf.gz",
        "combined_gvcfs/GRCh38.p12.genome.XXonly.chr17.gatk.combinegvcf.g.vcf.gz",
        "combined_gvcfs/GRCh38.p12.genome.XXonly.chr18.gatk.combinegvcf.g.vcf.gz",
        "combined_gvcfs/GRCh38.p12.genome.XXonly.chr19.gatk.combinegvcf.g.vcf.gz",
        "combined_gvcfs/GRCh38.p12.genome.XXonly.chr20.gatk.combinegvcf.g.vcf.gz",
        "combined_gvcfs/GRCh38.p12.genome.XXonly.chr21.gatk.combinegvcf.g.vcf.gz",
        "combined_gvcfs/GRCh38.p12.genome.XXonly.chr22.gatk.combinegvcf.g.vcf.gz",
        "combined_gvcfs/GRCh38.p12.genome.XXonly.chrX.gatk.combinegvcf.g.vcf.gz"
    input: #generate gvcf
        expand("gvcfs/{rna}/{rna}.GRCh38.p12.genome.XXonly.chr{chr}.g.vcf.gz", rna=config["all_rna_samples"], chr=chr_to_genotype)


# GATK on autosomes and X chromosomes
rule gatk_gvcf:
	input:
		ref = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa",
		bam = "/data/storage/SAYRES/placenta_YPOPS/06_HISAT_RNA/Aligned_BAMS/stranded_RF/{rna}_HISAT_pair_trim_sort_mkdup_rdgrp_XX.bam",
		bai = "/data/storage/SAYRES/placenta_YPOPS/06_HISAT_RNA/Aligned_BAMS/stranded_RF/{rna}_HISAT_pair_trim_sort_mkdup_rdgrp_XX.bam.bai"
	output:
		"gvcfs/{rna}/{rna}.GRCh38.p12.genome.XXonly.chr{chr}.g.vcf.gz"
	params:
		gatk = gatk_path,
		chrm_n = "chr{chr}"
	threads:
		4
	shell:
		"{params.gatk} "
		"HaplotypeCaller -R {input.ref} -I {input.bam} -L {params.chrm_n} "
		"--emit-ref-confidence GVCF --output {output}"

rule gatk_combinegvcfs_chr1:
	input:
		ref = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa",
		gvcfs = expand(
			"gvcfs/{rna}/{rna}.GRCh38.p12.genome.XXonly.chr1.g.vcf.gz", rna=config["all_rna_samples"])
	output:
		"combined_gvcfs/GRCh38.p12.genome.XXonly.chr1.gatk.combinegvcf.g.vcf.gz"
	params:
		gatk = gatk_path,
		chr_n = "chr1"
	threads:
		4
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell(
			"""{params.gatk} --java-options "-Xmx10g" """
			"""CombineGVCFs -R {input.ref} {variant_files} --intervals {params.chr_n} -O {output}""")

rule gatk_combinegvcfs_chr2:
	input:
		ref = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa",
		gvcfs = expand(
			"gvcfs/{rna}/{rna}.GRCh38.p12.genome.XXonly.chr2.g.vcf.gz", rna=config["all_rna_samples"])
	output:
		"combined_gvcfs/GRCh38.p12.genome.XXonly.chr2.gatk.combinegvcf.g.vcf.gz"
	params:
		gatk = gatk_path,
		chr_n = "chr2"
	threads:
		4
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell(
			"""{params.gatk} --java-options "-Xmx10g" """
			"""CombineGVCFs -R {input.ref} {variant_files} --intervals {params.chr_n} -O {output}""")

rule gatk_combinegvcfs_chr3:
	input:
		ref = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa",
		gvcfs = expand(
			"gvcfs/{rna}/{rna}.GRCh38.p12.genome.XXonly.chr3.g.vcf.gz", rna=config["all_rna_samples"])
	output:
		"combined_gvcfs/GRCh38.p12.genome.XXonly.chr3.gatk.combinegvcf.g.vcf.gz"
	params:
		gatk = gatk_path,
		chr_n = "chr3"
	threads:
		4
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell(
			"""{params.gatk} --java-options "-Xmx10g" """
			"""CombineGVCFs -R {input.ref} {variant_files} --intervals {params.chr_n} -O {output}""")

rule gatk_combinegvcfs_chr4:
	input:
		ref = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa",
		gvcfs = expand(
			"gvcfs/{rna}/{rna}.GRCh38.p12.genome.XXonly.chr4.g.vcf.gz", rna=config["all_rna_samples"])
	output:
		"combined_gvcfs/GRCh38.p12.genome.XXonly.chr4.gatk.combinegvcf.g.vcf.gz"
	params:
		gatk = gatk_path,
		chr_n = "chr4"
	threads:
		4
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell(
			"""{params.gatk} --java-options "-Xmx10g" """
			"""CombineGVCFs -R {input.ref} {variant_files} --intervals {params.chr_n} -O {output}""")

rule gatk_combinegvcfs_chr5:
	input:
		ref = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa",
		gvcfs = expand(
			"gvcfs/{rna}/{rna}.GRCh38.p12.genome.XXonly.chr5.g.vcf.gz", rna=config["all_rna_samples"])
	output:
		"combined_gvcfs/GRCh38.p12.genome.XXonly.chr5.gatk.combinegvcf.g.vcf.gz"
	params:
		gatk = gatk_path,
		chr_n = "chr5"
	threads:
		4
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell(
			"""{params.gatk} --java-options "-Xmx10g" """
			"""CombineGVCFs -R {input.ref} {variant_files} --intervals {params.chr_n} -O {output}""")

rule gatk_combinegvcfs_chr6:
	input:
		ref = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa",
		gvcfs = expand(
			"gvcfs/{rna}/{rna}.GRCh38.p12.genome.XXonly.chr6.g.vcf.gz", rna=config["all_rna_samples"])
	output:
		"combined_gvcfs/GRCh38.p12.genome.XXonly.chr6.gatk.combinegvcf.g.vcf.gz"
	params:
		gatk = gatk_path,
		chr_n = "chr6"
	threads:
		4
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell(
			"""{params.gatk} --java-options "-Xmx10g" """
			"""CombineGVCFs -R {input.ref} {variant_files} --intervals {params.chr_n} -O {output}""")

rule gatk_combinegvcfs_chr7:
	input:
		ref = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa",
		gvcfs = expand(
			"gvcfs/{rna}/{rna}.GRCh38.p12.genome.XXonly.chr7.g.vcf.gz", rna=config["all_rna_samples"])
	output:
		"combined_gvcfs/GRCh38.p12.genome.XXonly.chr7.gatk.combinegvcf.g.vcf.gz"
	params:
		gatk = gatk_path,
		chr_n = "chr7"
	threads:
		4
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell(
			"""{params.gatk} --java-options "-Xmx10g" """
			"""CombineGVCFs -R {input.ref} {variant_files} --intervals {params.chr_n} -O {output}""")

rule gatk_combinegvcfs_chr8:
	input:
		ref = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa",
		gvcfs = expand(
			"gvcfs/{rna}/{rna}.GRCh38.p12.genome.XXonly.chr8.g.vcf.gz", rna=config["all_rna_samples"])
	output:
		"combined_gvcfs/GRCh38.p12.genome.XXonly.chr8.gatk.combinegvcf.g.vcf.gz"
	params:
		gatk = gatk_path,
		chr_n = "chr8"
	threads:
		4
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell(
			"""{params.gatk} --java-options "-Xmx10g" """
			"""CombineGVCFs -R {input.ref} {variant_files} --intervals {params.chr_n} -O {output}""")

rule gatk_combinegvcfs_chr9:
	input:
		ref = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa",
		gvcfs = expand(
			"gvcfs/{rna}/{rna}.GRCh38.p12.genome.XXonly.chr9.g.vcf.gz", rna=config["all_rna_samples"])
	output:
		"combined_gvcfs/GRCh38.p12.genome.XXonly.chr9.gatk.combinegvcf.g.vcf.gz"
	params:
		gatk = gatk_path,
		chr_n = "chr9"
	threads:
		4
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell(
			"""{params.gatk} --java-options "-Xmx10g" """
			"""CombineGVCFs -R {input.ref} {variant_files} --intervals {params.chr_n} -O {output}""")

rule gatk_combinegvcfs_chr10:
	input:
		ref = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa",
		gvcfs = expand(
			"gvcfs/{rna}/{rna}.GRCh38.p12.genome.XXonly.chr10.g.vcf.gz", rna=config["all_rna_samples"])
	output:
		"combined_gvcfs/GRCh38.p12.genome.XXonly.chr10.gatk.combinegvcf.g.vcf.gz"
	params:
		gatk = gatk_path,
		chr_n = "chr10"
	threads:
		4
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell(
			"""{params.gatk} --java-options "-Xmx10g" """
			"""CombineGVCFs -R {input.ref} {variant_files} --intervals {params.chr_n} -O {output}""")

rule gatk_combinegvcfs_chr11:
	input:
		ref = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa",
		gvcfs = expand(
			"gvcfs/{rna}/{rna}.GRCh38.p12.genome.XXonly.chr11.g.vcf.gz", rna=config["all_rna_samples"])
	output:
		"combined_gvcfs/GRCh38.p12.genome.XXonly.chr11.gatk.combinegvcf.g.vcf.gz"
	params:
		gatk = gatk_path,
		chr_n = "chr11"
	threads:
		4
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell(
			"""{params.gatk} --java-options "-Xmx10g" """
			"""CombineGVCFs -R {input.ref} {variant_files} --intervals {params.chr_n} -O {output}""")

rule gatk_combinegvcfs_chr12:
	input:
		ref = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa",
		gvcfs = expand(
			"gvcfs/{rna}/{rna}.GRCh38.p12.genome.XXonly.chr12.g.vcf.gz", rna=config["all_rna_samples"])
	output:
		"combined_gvcfs/GRCh38.p12.genome.XXonly.chr12.gatk.combinegvcf.g.vcf.gz"
	params:
		gatk = gatk_path,
		chr_n = "chr12"
	threads:
		4
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell(
			"""{params.gatk} --java-options "-Xmx10g" """
			"""CombineGVCFs -R {input.ref} {variant_files} --intervals {params.chr_n} -O {output}""")

rule gatk_combinegvcfs_chr13:
	input:
		ref = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa",
		gvcfs = expand(
			"gvcfs/{rna}/{rna}.GRCh38.p12.genome.XXonly.chr13.g.vcf.gz", rna=config["all_rna_samples"])
	output:
		"combined_gvcfs/GRCh38.p12.genome.XXonly.chr13.gatk.combinegvcf.g.vcf.gz"
	params:
		gatk = gatk_path,
		chr_n = "chr13"
	threads:
		4
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell(
			"""{params.gatk} --java-options "-Xmx10g" """
			"""CombineGVCFs -R {input.ref} {variant_files} --intervals {params.chr_n} -O {output}""")

rule gatk_combinegvcfs_chr14:
	input:
		ref = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa",
		gvcfs = expand(
			"gvcfs/{rna}/{rna}.GRCh38.p12.genome.XXonly.chr14.g.vcf.gz", rna=config["all_rna_samples"])
	output:
		"combined_gvcfs/GRCh38.p12.genome.XXonly.chr14.gatk.combinegvcf.g.vcf.gz"
	params:
		gatk = gatk_path,
		chr_n = "chr14"
	threads:
		4
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell(
			"""{params.gatk} --java-options "-Xmx10g" """
			"""CombineGVCFs -R {input.ref} {variant_files} --intervals {params.chr_n} -O {output}""")

rule gatk_combinegvcfs_chr15:
	input:
		ref = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa",
		gvcfs = expand(
			"gvcfs/{rna}/{rna}.GRCh38.p12.genome.XXonly.chr15.g.vcf.gz", rna=config["all_rna_samples"])
	output:
		"combined_gvcfs/GRCh38.p12.genome.XXonly.chr15.gatk.combinegvcf.g.vcf.gz"
	params:
		gatk = gatk_path,
		chr_n = "chr15"
	threads:
		4
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell(
			"""{params.gatk} --java-options "-Xmx10g" """
			"""CombineGVCFs -R {input.ref} {variant_files} --intervals {params.chr_n} -O {output}""")

rule gatk_combinegvcfs_chr16:
	input:
		ref = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa",
		gvcfs = expand(
			"gvcfs/{rna}/{rna}.GRCh38.p12.genome.XXonly.chr16.g.vcf.gz", rna=config["all_rna_samples"])
	output:
		"combined_gvcfs/GRCh38.p12.genome.XXonly.chr16.gatk.combinegvcf.g.vcf.gz"
	params:
		gatk = gatk_path,
		chr_n = "chr16"
	threads:
		4
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell(
			"""{params.gatk} --java-options "-Xmx10g" """
			"""CombineGVCFs -R {input.ref} {variant_files} --intervals {params.chr_n} -O {output}""")

rule gatk_combinegvcfs_chr17:
	input:
		ref = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa",
		gvcfs = expand(
			"gvcfs/{rna}/{rna}.GRCh38.p12.genome.XXonly.chr17.g.vcf.gz", rna=config["all_rna_samples"])
	output:
		"combined_gvcfs/GRCh38.p12.genome.XXonly.chr17.gatk.combinegvcf.g.vcf.gz"
	params:
		gatk = gatk_path,
		chr_n = "chr17"
	threads:
		4
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell(
			"""{params.gatk} --java-options "-Xmx10g" """
			"""CombineGVCFs -R {input.ref} {variant_files} --intervals {params.chr_n} -O {output}""")

rule gatk_combinegvcfs_chr18:
	input:
		ref = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa",
		gvcfs = expand(
			"gvcfs/{rna}/{rna}.GRCh38.p12.genome.XXonly.chr18.g.vcf.gz", rna=config["all_rna_samples"])
	output:
		"combined_gvcfs/GRCh38.p12.genome.XXonly.chr18.gatk.combinegvcf.g.vcf.gz"
	params:
		gatk = gatk_path,
		chr_n = "chr18"
	threads:
		4
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell(
			"""{params.gatk} --java-options "-Xmx10g" """
			"""CombineGVCFs -R {input.ref} {variant_files} --intervals {params.chr_n} -O {output}""")

rule gatk_combinegvcfs_chr19:
	input:
		ref = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa",
		gvcfs = expand(
			"gvcfs/{rna}/{rna}.GRCh38.p12.genome.XXonly.chr19.g.vcf.gz", rna=config["all_rna_samples"])
	output:
		"combined_gvcfs/GRCh38.p12.genome.XXonly.chr19.gatk.combinegvcf.g.vcf.gz"
	params:
		gatk = gatk_path,
		chr_n = "chr19"
	threads:
		4
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell(
			"""{params.gatk} --java-options "-Xmx10g" """
			"""CombineGVCFs -R {input.ref} {variant_files} --intervals {params.chr_n} -O {output}""")

rule gatk_combinegvcfs_chr20:
	input:
		ref = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa",
		gvcfs = expand(
			"gvcfs/{rna}/{rna}.GRCh38.p12.genome.XXonly.chr20.g.vcf.gz", rna=config["all_rna_samples"])
	output:
		"combined_gvcfs/GRCh38.p12.genome.XXonly.chr20.gatk.combinegvcf.g.vcf.gz"
	params:
		gatk = gatk_path,
		chr_n = "chr20"
	threads:
		4
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell(
			"""{params.gatk} --java-options "-Xmx10g" """
			"""CombineGVCFs -R {input.ref} {variant_files} --intervals {params.chr_n} -O {output}""")

rule gatk_combinegvcfs_chr21:
	input:
		ref = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa",
		gvcfs = expand(
			"gvcfs/{rna}/{rna}.GRCh38.p12.genome.XXonly.chr21.g.vcf.gz", rna=config["all_rna_samples"])
	output:
		"combined_gvcfs/GRCh38.p12.genome.XXonly.chr21.gatk.combinegvcf.g.vcf.gz"
	params:
		gatk = gatk_path,
		chr_n = "chr21"
	threads:
		4
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell(
			"""{params.gatk} --java-options "-Xmx10g" """
			"""CombineGVCFs -R {input.ref} {variant_files} --intervals {params.chr_n} -O {output}""")

rule gatk_combinegvcfs_chr22:
	input:
		ref = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa",
		gvcfs = expand(
			"gvcfs/{rna}/{rna}.GRCh38.p12.genome.XXonly.chr22.g.vcf.gz", rna=config["all_rna_samples"])
	output:
		"combined_gvcfs/GRCh38.p12.genome.XXonly.chr22.gatk.combinegvcf.g.vcf.gz"
	params:
		gatk = gatk_path,
		chr_n = "chr22"
	threads:
		4
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell(
			"""{params.gatk} --java-options "-Xmx10g" """
			"""CombineGVCFs -R {input.ref} {variant_files} --intervals {params.chr_n} -O {output}""")

rule gatk_combinegvcfs_chrX:
	input:
		ref = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa",
		gvcfs = expand(
			"gvcfs/{rna}/{rna}.GRCh38.p12.genome.XXonly.chrX.g.vcf.gz", rna=config["all_rna_samples"])
	output:
		"combined_gvcfs/GRCh38.p12.genome.XXonly.chrX.gatk.combinegvcf.g.vcf.gz"
	params:
		gatk = gatk_path,
		chr_n = "chrX"
	threads:
		4
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell(
			"""{params.gatk} --java-options "-Xmx10g" """
			"""CombineGVCFs -R {input.ref} {variant_files} --intervals {params.chr_n} -O {output}""")

rule gatk_genotypegvcf:
	input:
		ref = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa",
		gvcf = "combined_gvcfs/GRCh38.p12.genome.XXonly.chr{chr}.gatk.combinegvcf.g.vcf.gz"
	output:
		"genotyped_vcfs/GRCh38.p12.genome.XXonly.chr{chr}.gatk.called.raw.vcf.gz"
	params:
		gatk = gatk_path
	threads:
		4
	shell:
		"""{params.gatk} --java-options "-Xmx10g" """
		"""GenotypeGVCFs -R {input.ref} -V {input.gvcf} -O {output}"""

#----------------------------------------
# Further processing VCF prior to VQSR
# 1. Restrict to biallelic sites
# 2. Subset VCF files for each individual
# 3. Keep only the heterozygous sites
# Do this for chr8 and chrX
#----------------------------------------
rule gatk_selectbiallelic_before_vqsr:
    input:
        ref = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa",
        vcf = "genotyped_vcfs/GRCh38.p12.genome.XXonly.chr{chr}.gatk.called.raw.vcf.gz"
    output:
        "genotyped_vcfs/chr{chr}.gatk.called.raw.biallelic.snp.vcf.gz"
    params:
        gatk = gatk_path
    shell:
        """{params.gatk} SelectVariants """
        """-R {input.ref} """
        """-V {input.vcf} """
        """-O {output} """
        """--select-type-to-include SNP """
        """--restrict-alleles-to BIALLELIC """

# rule subset_individuals_before_vqsr:
#     input:
#         "genotyped_vcfs/chr{chr}.gatk.called.raw.biallelic.snp.vcf.gz"
#     output:
#         "genotyped_vcfs/chr{chr}.gatk.called.raw.biallelic.snp.{rna}.vcf"
#     params:
#         bcftools = bcftools_path,
#         rna = "{rna}"
#     shell:
#         """{params.bcftools} view -s {params.rna} {input} > {output}"""

# After subsetting for each individual. In some individuals,
# the genotypes could be homozygous for the reference. This next rule is to remove these sites.
rule gatk_selectheterozygous_before_vqsr:
    input:
        ref = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa",
        vcf = "genotyped_vcfs/chr{chr}.gatk.called.raw.biallelic.snp.{rna}.vcf"
    output:
        "genotyped_vcfs/chr{chr}.gatk.called.raw.biallelic.snp.{rna}.het.vcf"
    params:
        gatk = gatk_path
    shell:
        """{params.gatk} SelectVariants """
        """-R {input.ref} """
        """-V {input.vcf} """
        """-O {output} """
        """-select "AC == 1" """
