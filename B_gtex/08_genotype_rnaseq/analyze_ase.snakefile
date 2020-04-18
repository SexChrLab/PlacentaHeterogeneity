import os

configfile: "asereadcounter_config.json"

chromosomes = ["8", "X"]

tissues = ["AdrenalGland", "ArteryAorta", "BrainAmygdala", "BrainCortex", "BreastMammaryTissue", "ColonSigmoid", "HeartLeftVentricle", "KidneyCortex", "Liver", "Lung", "MuscleSkeletal", "NerveTibial", "Ovary", "Pancreas", "Pituitary", "Spleen", "Stomach", "Thyroid", "Uterus", "Vagina"]

combiC = []
for key in config["dna_rna"]:
    for item in config["dna_rna"][key]:
        combiC.append((key, item))

import itertools
combiList=list()
for c in combiC:
    combiList.append(c[0]+"_"+c[1])

rule all:
    input: #calc_allele_balance
        expand("asereadcounter/HISAT/after_vqsr/chr{chr}/{combo}_chr{chr}_allele_balance.tsv", combo=combiList, chr=chromosomes),
        expand("asereadcounter/HISAT/after_vqsr/chr8/{tissue}_chr8_median_allele_balance.txt", tissue=tissues), #chr8
        expand("asereadcounter/HISAT/after_vqsr/chrX/{tissue}_chrX_nonpars_median_allele_balance.txt", tissue=tissues)

rule calc_allele_balance:
    input:
        "asereadcounter/HISAT/after_vqsr/chr{chr}/{dna}_{rna}_chr{chr}.tsv"
    output:
        "asereadcounter/HISAT/after_vqsr/chr{chr}/{dna}_{rna}_chr{chr}_allele_balance.tsv"
    params:
        script = "/scratch/tphung3/PlacentaSexDiff/B_gtex/07_analyze_ase/scripts/calc_allele_balance.py"
    shell:
        """
        python {params.script} --input {input} --output {output}
        """

rule calc_median_allele_balance_per_tissue_chr8:
    output:
        "asereadcounter/HISAT/after_vqsr/chr8/{tissue}_chr8_median_allele_balance.txt"
    params:
        script = "/scratch/tphung3/PlacentaSexDiff/B_gtex/07_analyze_ase/scripts/calc_median_allele_balance_per_tissue.py",
        tissue = "{tissue}",
        chr = "chr8",
        dir = "asereadcounter/HISAT/after_vqsr/chr8",
        threshold = 10
    shell:
        """
        python {params.script} --tissue {params.tissue} --chromosome {params.chr} --files_dir {params.dir} --threshold {params.threshold} --include_pars yes --out {output}
        """

rule calc_median_allele_balance_per_tissue_chrX_nonpars:
    output:
        "asereadcounter/HISAT/after_vqsr/chrX/{tissue}_chrX_nonpars_median_allele_balance.txt"
    params:
        script = "/scratch/tphung3/PlacentaSexDiff/B_gtex/07_analyze_ase/scripts/calc_median_allele_balance_per_tissue.py",
        tissue = "{tissue}",
        chr = "chrX",
        dir = "asereadcounter/HISAT/after_vqsr/chrX",
        threshold = 10
    shell:
        """
        python {params.script} --tissue {params.tissue} --chromosome {params.chr} --files_dir {params.dir} --threshold {params.threshold} --include_pars no --out {output}
        """
