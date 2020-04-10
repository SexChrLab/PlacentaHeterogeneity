import os

configfile: "process_dna_females_config.json"

adapter_path = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/adapter_sequence.fa"

chr_to_genotype = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"]
chr_current = ["8", "X"]

gatk_path = "/home/tphung3/softwares/gatk-4.1.0.0/gatk"
bcftools_path = "bcftools"


rule all:
    # Further processing
    input:
        expand("vqsr/chr{chr}.gatk.called.vqsr.sv.biallelic.snp.{sample}.het.vcf", chr=chr_current, sample=config["all_dna_samples"])
    input:
        expand("vqsr/chr{chr}.gatk.called.vqsr.sv.biallelic.snp.{sample}.vcf", chr=chr_current, sample=config["all_dna_samples"])
    input: #vqsr
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
        expand("gvcfs/{sample}/{sample}.GRCh38.p12.genome.XXonly.chr{chr}.g.vcf.gz", sample=config["all_dna_samples"], chr=chr_to_genotype)
    input: #trim, map, index, and bam stats
        expand("processed_bams/dna/{sample}.GRCh38.p12.genome.XXonly.mkdup.sorted.mkdup.bam.bai", sample=config["all_dna_samples"]),
        expand("bam_stats/dna/{sample}.GRCh38.p12.genome.XXonly.mkdup.sorted.stats", sample=config["all_dna_samples"]),
        expand("picard_stats/{sample}.GRCh38.p12.genome.XXonly.picard_mkdup_metrics.txt", sample=config["all_dna_samples"])


rule trim_adapters_paired_bbduk_dna:
    input:
        fq_1 = lambda wildcards: os.path.join(config[wildcards.sample]["fq_path"], config[wildcards.sample]["fq_1"]),
        fq_2 = lambda wildcards: os.path.join(config[wildcards.sample]["fq_path"], config[wildcards.sample]["fq_2"])
    output:
        out_fq_1 = "trimmed_fastqs_dna/{sample}_trimmed_read1.fastq.gz",
        out_fq_2 = "trimmed_fastqs_dna/{sample}_trimmed_read2.fastq.gz"
    params:
        adapter = adapter_path
    threads:
        2
    shell:
        "bbduk.sh -Xmx3g in1={input.fq_1} in2={input.fq_2} out1={output.out_fq_1} out2={output.out_fq_2} ref={params.adapter} qtrim=rl trimq=30 minlen=75 maq=20"

rule map_and_process_trimmed_reads_dna:
	input:
		fq_1 = "trimmed_fastqs_dna/{sample}_trimmed_read1.fastq.gz",
		fq_2 = "trimmed_fastqs_dna/{sample}_trimmed_read2.fastq.gz",
		fai_xx = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa.fai",
		ref_xx = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa"
	output:
		"processed_bams/dna/{sample}.GRCh38.p12.genome.XXonly.mkdup.sorted.bam"
	params:
		id = lambda wildcards: config[wildcards.sample]["ID"],
		sm = lambda wildcards: config[wildcards.sample]["SM"],
		lb = lambda wildcards: config[wildcards.sample]["LB"],
		pu = lambda wildcards: config[wildcards.sample]["PU"],
		pl = lambda wildcards: config[wildcards.sample]["PL"]
	threads:
		4
	shell:
		"bwa mem -t {threads} -R "
		"'@RG\\tID:{params.id}\\tSM:{params.sm}\\tLB:{params.lb}\\tPU:{params.pu}\\tPL:{params.pl}' "
		"{input.ref_xx} {input.fq_1} {input.fq_2} "
		"| samblaster "
		"| samtools fixmate -O bam - - | samtools sort "
		"-O bam -o {output}"

rule index_bam_dna:
	input:
		"processed_bams/dna/{sample}.GRCh38.p12.genome.XXonly.mkdup.sorted.bam"
	output:
		"processed_bams/dna/{sample}.GRCh38.p12.genome.XXonly.mkdup.sorted.bam.bai"
	shell:
		"samtools index {input}"

rule bam_stats_dna:
	input:
		"processed_bams/dna/{sample}.GRCh38.p12.genome.XXonly.mkdup.sorted.bam"
	output:
		"bam_stats/dna/{sample}.GRCh38.p12.genome.XXonly.mkdup.sorted.stats"
	shell:
		"samtools stats {input} | grep ^SN | cut -f 2- > {output}"

rule picard_mkdups:
    input:
        bam = "processed_bams/dna/{sample}.GRCh38.p12.genome.XXonly.mkdup.sorted.bam",
        bai = "processed_bams/dna/{sample}.GRCh38.p12.genome.XXonly.mkdup.sorted.bam.bai"
    output:
        bam = "processed_bams/dna/{sample}.GRCh38.p12.genome.XXonly.mkdup.sorted.mkdup.bam",
        metrics = "picard_stats/{sample}.GRCh38.p12.genome.XXonly.picard_mkdup_metrics.txt"
    threads: 4
    shell:
        "picard -Xmx14g MarkDuplicates I={input.bam} O={output.bam} "
        "M={output.metrics}"

rule index_mkdup_bams_dna:
    input:
        "processed_bams/dna/{sample}.GRCh38.p12.genome.XXonly.mkdup.sorted.mkdup.bam"
    output:
        "processed_bams/dna/{sample}.GRCh38.p12.genome.XXonly.mkdup.sorted.mkdup.bam.bai"
    shell:
        "samtools index {input}"

# GATK on autosomes and X chromosomes
rule gatk_gvcf:
	input:
		ref = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa",
		bam = "processed_bams/dna/{sample}.GRCh38.p12.genome.XXonly.mkdup.sorted.mkdup.bam",
		bai = "processed_bams/dna/{sample}.GRCh38.p12.genome.XXonly.mkdup.sorted.mkdup.bam.bai"
	output:
		"gvcfs/{sample}/{sample}.GRCh38.p12.genome.XXonly.chr{chr}.g.vcf.gz"
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
			"gvcfs/{sample}/{sample}.GRCh38.p12.genome.XXonly.chr1.g.vcf.gz", sample=config["all_dna_samples"])
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
			"gvcfs/{sample}/{sample}.GRCh38.p12.genome.XXonly.chr2.g.vcf.gz", sample=config["all_dna_samples"])
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
			"gvcfs/{sample}/{sample}.GRCh38.p12.genome.XXonly.chr3.g.vcf.gz", sample=config["all_dna_samples"])
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
			"gvcfs/{sample}/{sample}.GRCh38.p12.genome.XXonly.chr4.g.vcf.gz", sample=config["all_dna_samples"])
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
			"gvcfs/{sample}/{sample}.GRCh38.p12.genome.XXonly.chr5.g.vcf.gz", sample=config["all_dna_samples"])
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
			"gvcfs/{sample}/{sample}.GRCh38.p12.genome.XXonly.chr6.g.vcf.gz", sample=config["all_dna_samples"])
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
			"gvcfs/{sample}/{sample}.GRCh38.p12.genome.XXonly.chr7.g.vcf.gz", sample=config["all_dna_samples"])
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
			"gvcfs/{sample}/{sample}.GRCh38.p12.genome.XXonly.chr8.g.vcf.gz", sample=config["all_dna_samples"])
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
			"gvcfs/{sample}/{sample}.GRCh38.p12.genome.XXonly.chr9.g.vcf.gz", sample=config["all_dna_samples"])
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
			"gvcfs/{sample}/{sample}.GRCh38.p12.genome.XXonly.chr10.g.vcf.gz", sample=config["all_dna_samples"])
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
			"gvcfs/{sample}/{sample}.GRCh38.p12.genome.XXonly.chr11.g.vcf.gz", sample=config["all_dna_samples"])
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
			"gvcfs/{sample}/{sample}.GRCh38.p12.genome.XXonly.chr12.g.vcf.gz", sample=config["all_dna_samples"])
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
			"gvcfs/{sample}/{sample}.GRCh38.p12.genome.XXonly.chr13.g.vcf.gz", sample=config["all_dna_samples"])
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
			"gvcfs/{sample}/{sample}.GRCh38.p12.genome.XXonly.chr14.g.vcf.gz", sample=config["all_dna_samples"])
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
			"gvcfs/{sample}/{sample}.GRCh38.p12.genome.XXonly.chr15.g.vcf.gz", sample=config["all_dna_samples"])
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
			"gvcfs/{sample}/{sample}.GRCh38.p12.genome.XXonly.chr16.g.vcf.gz", sample=config["all_dna_samples"])
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
			"gvcfs/{sample}/{sample}.GRCh38.p12.genome.XXonly.chr17.g.vcf.gz", sample=config["all_dna_samples"])
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
			"gvcfs/{sample}/{sample}.GRCh38.p12.genome.XXonly.chr18.g.vcf.gz", sample=config["all_dna_samples"])
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
			"gvcfs/{sample}/{sample}.GRCh38.p12.genome.XXonly.chr19.g.vcf.gz", sample=config["all_dna_samples"])
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
			"gvcfs/{sample}/{sample}.GRCh38.p12.genome.XXonly.chr20.g.vcf.gz", sample=config["all_dna_samples"])
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
			"gvcfs/{sample}/{sample}.GRCh38.p12.genome.XXonly.chr21.g.vcf.gz", sample=config["all_dna_samples"])
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
			"gvcfs/{sample}/{sample}.GRCh38.p12.genome.XXonly.chr22.g.vcf.gz", sample=config["all_dna_samples"])
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
			"gvcfs/{sample}/{sample}.GRCh38.p12.genome.XXonly.chrX.g.vcf.gz", sample=config["all_dna_samples"])
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
        """-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an InbreedingCoeff """
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
        """-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an InbreedingCoeff """
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
#         "vqsr/chr{chr}.gatk.called.vqsr.sv.biallelic.snp.{sample}.vcf"
#     params:
#         bcftools = bcftools_path,
#         sample = "{sample}"
#     shell:
#         """{params.bcftools} view -s {params.sample} {input} > {output}"""

# After subsetting for each individual. In some individuals,
# the genotypes could be homozygous for the reference. This next rule is to remove these sites.
rule gatk_selectheterozygous:
    input:
        ref = "/scratch/tphung3/PlacentaSexDiff/A_placenta/01_process_dna/refs/GRCh38.p12.genome.XXonly.fa",
        vcf = "vqsr/chr{chr}.gatk.called.vqsr.sv.biallelic.snp.{sample}.vcf"
    output:
        "vqsr/chr{chr}.gatk.called.vqsr.sv.biallelic.snp.{sample}.het.vcf"
    params:
        gatk = gatk_path
    shell:
        """{params.gatk} SelectVariants """
        """-R {input.ref} """
        """-V {input.vcf} """
        """-O {output} """
        """-select "AC == 1" """
