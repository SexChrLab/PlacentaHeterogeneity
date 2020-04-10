import os

# This snakefile extracts read groups from the bam file
configfile: "/scratch/tphung3/PlacentaSexDiff/B_gtex/04_process_rna/process_rna_config.json"

rule all:
    input:
        expand("read_groups/rna/{sample}_RG.txt", sample=config["rna_samples"])

rule extract_readgroups_rna:
    input:
        BAM = os.path.join(config["GTEX_bam_dir"], "{sample}.Aligned.sortedByCoord.out.patched.md.bam")
    output:
        "read_groups/rna/{sample}_RG.txt"
    shell:
        """
        samtools view -H {input.BAM} | grep "^@RG" > {output}
        """
