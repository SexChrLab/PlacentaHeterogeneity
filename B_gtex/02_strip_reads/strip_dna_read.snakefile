import os

configfile: "strip_dna_read_config.json"

rule all:
    input:
        expand("strip_reads/{sample}/{sample}_xyalign.log", sample=config["all_dna_samples"])

rule strip_reads:
    input:
        BAM = os.path.join(config["GTEX_bam_dir"], "{sample}.bam")
    output:
        LOG = "strip_reads/{sample}/{sample}_xyalign.log"
    conda:
        "/scratch/tphung3/PlacentaSexDiff/B_gtex/02_strip_reads/envs/xyalign.yml"
    params:
        DIR = "strip_reads/{sample}",
        SAMPLE_ID = "{sample}",
        cpus="4",
        xmx="4g",
        compression="3"
    shell:
        "xyalign --STRIP_READS --ref null --bam {input.BAM} --cpus {params.cpus} --xmx {params.xmx} --sample_id {params.SAMPLE_ID} --output_dir {params.DIR} --chromosomes ALL"
