import os

# This snakefile extracts read groups from the bam file
configfile: "/scratch/tphung3/PlacentaSexDiff/B_gtex/05_process_dna/process_dna_config.json"

rule all:
    input:
        expand("read_groups/dna/{sample}_RG.txt", sample=config["dna_samples"])

rule extract_readgroups_dna:
    input:
        BAM = os.path.join(config["GTEX_bam_dir"], "{sample}.bam")
    output:
        "read_groups/dna/{sample}_RG.txt"
    shell:
        """
        samtools view -H {input.BAM} | grep "^@RG" > {output}
        """
