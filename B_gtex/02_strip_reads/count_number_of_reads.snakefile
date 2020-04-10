import os

# configfile: "strip_rna_read_config.json"
configfile: "strip_dna_read_config.json"
rule all:
    input:
        expand("count_number_of_reads/{sample}_number_of_reads.txt", sample=config["dna_samples"])
rule count_number_of_reads:
    input:
    output:
        "count_number_of_reads/{sample}_number_of_reads.txt"
    params:
        id = "{sample}",
        script = "count_number_of_reads.py"
    shell:
        """
        python {params.script} {params.id} {output}
        """
