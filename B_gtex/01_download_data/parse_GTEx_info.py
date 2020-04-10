import pandas as pd
import os

# Load the file that has the female GTEx id
female_ids = set()
with open("c://Users/tuyen/Documents/postdoc_asu/projects/GTEx/GTEx_v8/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS_females.txt", "r") as f:
    for line in f:
        if not line.startswith("SUBJID"):
            female_ids.add(line.rstrip("\n").split("\t")[0])

attributes = pd.read_csv("c://Users/tuyen/Documents/postdoc_asu/projects/GTEx/GTEx_v8/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", sep="\t")

tissues = ['Adrenal Gland', 'Artery - Aorta', 'Bladder', 'Brain - Amygdala', 'Brain - Cortex', 'Breast - Mammary Tissue', 'Colon - Sigmoid', 'Heart - Left Ventricle', 'Kidney - Cortex', 'Liver', 'Lung', 'Muscle - Skeletal', 'Nerve - Tibial', 'Ovary', 'Pancreas', 'Pituitary', 'Spleen', 'Stomach', 'Thyroid', 'Uterus', 'Vagina']
# Subset for each tissue
for tissue in tissues:
    tissue_subset = attributes[attributes['SMTSD'] == tissue]

    # Save to a file
    outfile_path = os.path.join("c://Users/tuyen/Documents/postdoc_asu/projects/PlacentaSexDiff/B_gtex/01_download_data/tissues_subset/", "GTEx_Analysis_v8_Annotations_SampleAttributesDS_" + tissue + ".txt")
    tissue_subset.to_csv(outfile_path, sep="\t", index=False)

for tissue in tissues:
    outfile = open(os.path.join("c://Users/tuyen/Documents/postdoc_asu/projects/PlacentaSexDiff/B_gtex/01_download_data/tissues_subset/", "GTEx_Analysis_v8_Annotations_SampleAttributesDS_" + tissue + "_females.txt"), "w")
    with open(os.path.join("c://Users/tuyen/Documents/postdoc_asu/projects/PlacentaSexDiff/B_gtex/01_download_data/tissues_subset/", "GTEx_Analysis_v8_Annotations_SampleAttributesDS_" + tissue + ".txt"), "r") as f:
        for line in f:
            if line.startswith('SAMPID'):
                print (line, file=outfile)
            else:
                items = line.rstrip('\n').split('\t')
                i = items[0].split('-')
                id = i[0] + '-' + i[1]
                if id in female_ids:
                    print (line, file=outfile)
