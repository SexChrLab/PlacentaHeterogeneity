#!/bin/bash

for tissue in AdrenalGland ArteryAorta BrainAmygdala BrainCortex BreastMammaryTissue ColonSigmoid HeartLeftVentricle KidneyCortex Liver Lung MuscleSkeletal NerveTibial Ovary Pancreas Pituitary Spleen Stomach Thyroid Uterus Vagina
do

Rscript compare_wes_rnaseq.R ${tissue}

done
