import os

configfile: "genotype_rnaseq_config.json"

chr_to_genotype = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"]
chr_current = ["8", "X"]

gatk_path = "/home/tphung3/softwares/gatk-4.1.0.0/gatk"
bcftools_path = "bcftools"


rule all:
    input:
        expand("vqsr/chr{chr}.gatk.called.vqsr.sv.biallelic.snp.{rna}.het.vcf", chr=chr_current, rna=config["all_rna_samples_short"])
    input:
        expand("vqsr/chr{chr}.gatk.called.vqsr.sv.biallelic.snp.{rna}.vcf", chr=chr_current, rna=config["all_rna_samples_short"])
    input: #filter with VQSR
        "vqsr/chr8.gatk.called.vqsr.sv.vcf.gz",
        "vqsr/chrX.gatk.called.vqsr.sv.vcf.gz"
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
		bam = "/scratch/tphung3/PlacentaSexDiff/B_gtex/04_process_rna/processed_bams/rna/{rna}.GRCh38.p12.genome.XXonly.sorted.merged.bam",
		bai = "/scratch/tphung3/PlacentaSexDiff/B_gtex/04_process_rna/processed_bams/rna/{rna}.GRCh38.p12.genome.XXonly.sorted.merged.bam.bai"
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

# ----------------
# Filter with VQSR
# ----------------
# chr8
rule gatk_variantrecalibrator_chr8:
    input:
        ref = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa",
        vcf = "genotyped_vcfs/GRCh38.p12.genome.XXonly.chr8.gatk.called.raw.vcf.gz",
        hapmap = "/scratch/tphung3/PopulationReferenceAlignment/analyses/compare_hardfilter_vqsr/scripts/hapmap_3.3.hg38.vcf.gz",
        omni = "/scratch/tphung3/PopulationReferenceAlignment/analyses/compare_hardfilter_vqsr/scripts/1000G_omni2.5.hg38.vcf.gz",
        thousandG = "/scratch/tphung3/PopulationReferenceAlignment/analyses/compare_hardfilter_vqsr/scripts/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
        dbsnp = "/scratch/tphung3/PopulationReferenceAlignment/analyses/compare_hardfilter_vqsr/scripts/dbsnp_138.hg38.vcf.gz"
    output:
        recal = "vqsr/chr8_output.recal",
        tranches = "vqsr/chr8_output.tranches"
    params:
        gatk = gatk_path
    shell:
        """{params.gatk} --java-options "-Xmx16g" VariantRecalibrator """
        """-R {input.ref} -V {input.vcf}  """
        """--resource:hapmap,known=false,training=true,truth=true,prior=15.0 {input.hapmap} """
        """--resource:omni,known=false,training=true,truth=false,prior=12.0 {input.omni} """
        """--resource:1000G,known=false,training=true,truth=false,prior=10.0 {input.thousandG} """
        """--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {input.dbsnp} """
        """-an QD -an ReadPosRankSum -an FS -an SOR -an InbreedingCoeff """
        """-mode SNP """
        """-O {output.recal} """
        """--tranches-file {output.tranches} """

rule gatk_applyvqsr_chr8:
    input:
        ref = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa",
        vcf = "genotyped_vcfs/GRCh38.p12.genome.XXonly.chr8.gatk.called.raw.vcf.gz",
        tranches = "vqsr/chr8_output.tranches",
        recal = "vqsr/chr8_output.recal"
    output:
        "vqsr/chr8.gatk.called.vqsr.vcf.gz"
    params:
        gatk = gatk_path
    shell:
        """{params.gatk} --java-options "-Xmx16g" ApplyVQSR """
        """-R {input.ref} """
        """-V {input.vcf} """
        """-O {output} """
        """--truth-sensitivity-filter-level 99.0 """
        """--tranches-file {input.tranches} """
        """--recal-file {input.recal} """
        """-mode SNP """

rule gatk_selectvariants_chr8:
    input:
        ref = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa",
        vcf = "vqsr/chr8.gatk.called.vqsr.vcf.gz"
    output:
        "vqsr/chr8.gatk.called.vqsr.sv.vcf.gz"
    params:
        gatk = gatk_path
    shell:
        """{params.gatk} --java-options "-Xmx16g" SelectVariants """
        """-R {input.ref} """
        """-V {input.vcf} """
        """--exclude-filtered """
        """-O {output} """

# chrX
rule gatk_variantrecalibrator_chrX:
    input:
        ref = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa",
        vcf = "genotyped_vcfs/GRCh38.p12.genome.XXonly.chrX.gatk.called.raw.vcf.gz",
        hapmap = "/scratch/tphung3/PopulationReferenceAlignment/analyses/compare_hardfilter_vqsr/scripts/hapmap_3.3.hg38.vcf.gz",
        omni = "/scratch/tphung3/PopulationReferenceAlignment/analyses/compare_hardfilter_vqsr/scripts/1000G_omni2.5.hg38.vcf.gz",
        thousandG = "/scratch/tphung3/PopulationReferenceAlignment/analyses/compare_hardfilter_vqsr/scripts/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
        dbsnp = "/scratch/tphung3/PopulationReferenceAlignment/analyses/compare_hardfilter_vqsr/scripts/dbsnp_138.hg38.vcf.gz"
    output:
        recal = "vqsr/chrX_output.recal",
        tranches = "vqsr/chrX_output.tranches"
    params:
        gatk = gatk_path
    shell:
        """{params.gatk} --java-options "-Xmx16g" VariantRecalibrator """
        """-R {input.ref} -V {input.vcf}  """
        """--resource:hapmap,known=false,training=true,truth=true,prior=15.0 {input.hapmap} """
        """--resource:omni,known=false,training=true,truth=false,prior=12.0 {input.omni} """
        """--resource:1000G,known=false,training=true,truth=false,prior=10.0 {input.thousandG} """
        """--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {input.dbsnp} """
        """-an QD -an ReadPosRankSum -an FS -an SOR -an InbreedingCoeff """
        """-mode SNP """
        """-O {output.recal} """
        """--tranches-file {output.tranches} """

rule gatk_applyvqsr_chrX:
    input:
        ref = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa",
        vcf = "genotyped_vcfs/GRCh38.p12.genome.XXonly.chrX.gatk.called.raw.vcf.gz",
        tranches = "vqsr/chrX_output.tranches",
        recal = "vqsr/chrX_output.recal"
    output:
        "vqsr/chrX.gatk.called.vqsr.vcf.gz"
    params:
        gatk = gatk_path
    shell:
        """{params.gatk} --java-options "-Xmx16g" ApplyVQSR """
        """-R {input.ref} """
        """-V {input.vcf} """
        """-O {output} """
        """--truth-sensitivity-filter-level 99.0 """
        """--tranches-file {input.tranches} """
        """--recal-file {input.recal} """
        """-mode SNP """

rule gatk_selectvariants_chrX:
    input:
        ref = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa",
        vcf = "vqsr/chrX.gatk.called.vqsr.vcf.gz"
    output:
        "vqsr/chrX.gatk.called.vqsr.sv.vcf.gz"
    params:
        gatk = gatk_path
    shell:
        """{params.gatk} --java-options "-Xmx16g" SelectVariants """
        """-R {input.ref} """
        """-V {input.vcf} """
        """--exclude-filtered """
        """-O {output} """

#----------------------------------------
# Further processing VCF
# 1. Restrict to biallelic sites
# 2. Subset VCF files for each individual
# 3. Keep only the heterozygous sites
# Do this for chr8 and chrX
#----------------------------------------
rule gatk_selectbiallelic:
    input:
        ref = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa",
        vcf = "vqsr/chr{chr}.gatk.called.vqsr.sv.vcf.gz"
    output:
        "vqsr/chr{chr}.gatk.called.vqsr.sv.biallelic.snp.vcf.gz"
    params:
        gatk = gatk_path
    shell:
        """{params.gatk} SelectVariants """
        """-R {input.ref} """
        """-V {input.vcf} """
        """-O {output} """
        """--select-type-to-include SNP """
        """--restrict-alleles-to BIALLELIC """

# rule subset_individuals_vqsr:
#     input:
#         "vqsr/chr{chr}.gatk.called.vqsr.sv.biallelic.snp.vcf.gz"
#     output:
#         "vqsr/chr{chr}.gatk.called.vqsr.sv.biallelic.snp.{rna}.vcf"
#     params:
#         bcftools = bcftools_path,
#         sample = "{rna}"
#     shell:
#         """{params.bcftools} view -s {params.sample} {input} > {output}"""

# After subsetting for each individual. In some individuals,
# the genotypes could be homozygous for the reference. This next rule is to remove these sites.
rule gatk_selectheterozygous:
    input:
        ref = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa",
        vcf = "vqsr/chr{chr}.gatk.called.vqsr.sv.biallelic.snp.{rna}.vcf"
    output:
        "vqsr/chr{chr}.gatk.called.vqsr.sv.biallelic.snp.{rna}.het.vcf"
    params:
        gatk = gatk_path
    shell:
        """{params.gatk} SelectVariants """
        """-R {input.ref} """
        """-V {input.vcf} """
        """-O {output} """
        """-select "AC == 1" """
