library(ggplot2)
library(ggpubr)

# chrX
chrX_adrenal = read.csv('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chrX/AdrenalGland_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_artery = read.csv('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chrX/ArteryAorta_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_amygdala = read.csv('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chrX/BrainAmygdala_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_cortex = read.csv('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chrX/BrainCortex_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_breast = read.csv('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chrX/BreastMammaryTissue_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_colon = read.csv('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chrX/ColonSigmoid_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_heart = read.csv('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chrX/HeartLeftVentricle_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_kidney = read.csv('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chrX/KidneyCortex_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_liver = read.csv('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chrX/Liver_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_lung = read.csv('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chrX/Lung_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_muscle = read.csv('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chrX/MuscleSkeletal_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_nerve = read.csv('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chrX/NerveTibial_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_ovary = read.csv('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chrX/Ovary_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_pancreas = read.csv('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chrX/Pancreas_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_pituitary = read.csv('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chrX/Pituitary_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_spleen = read.csv('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chrX/Spleen_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_stomach = read.csv('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chrX/Stomach_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_thyroid = read.csv('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chrX/Thyroid_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_uterus = read.csv('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chrX/Uterus_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_vagina = read.csv('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chrX/Vagina_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_placenta = read.csv('/scratch/tphung3/PlacentaSexDiff/A_placenta/05_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chrX/placenta_chrX_median_allele_balance.txt', sep = '\t')

adrenal_lab = rep('Adrenal', 10)
artery_lab = rep('Artery', 10)
amygdala_lab = rep('Amygdala', 10)
cortex_lab = rep('Cortex', 10)
breast_lab = rep('Breast', 10)
colon_lab = rep('Colon', 10)
heart_lab = rep('Heart', 10)
kidney_lab = rep('Kidney', 10)
liver_lab = rep('Liver', 10)
lung_lab = rep('Lung', 10)
muscle_lab = rep('Muscle', 10)
nerve_lab = rep('Nerve', 10)
ovary_lab = rep('Ovary', 10)
pancreas_lab = rep('Pancreas', 10)
pituitary_lab = rep('Pituitary', 10)
spleen_lab = rep('Spleen', 10)
stomach_lab = rep('Stomach', 10)
thyroid_lab = rep('Thyroid', 10)
uterus_lab = rep('Uterus', 10)
vagina_lab = rep('Vagina', 10)
placenta_lab = rep('Placenta', 60)

gtex_lab = rep('GTEX', 200)

chrX_data = data.frame(median = c(chrX_adrenal$median_allele_balance_subset, 
                                  chrX_artery$median_allele_balance_subset, 
                                  chrX_amygdala$median_allele_balance_subset, 
                                  chrX_cortex$median_allele_balance_subset, 
                                  chrX_breast$median_allele_balance_subset, 
                                  chrX_colon$median_allele_balance_subset, 
                                  chrX_heart$median_allele_balance_subset, 
                                  chrX_kidney$median_allele_balance_subset, 
                                  chrX_liver$median_allele_balance_subset, 
                                  chrX_lung$median_allele_balance_subset, 
                                  chrX_muscle$median_allele_balance_subset, 
                                  chrX_nerve$median_allele_balance_subset, 
                                  chrX_ovary$median_allele_balance_subset, 
                                  chrX_pancreas$median_allele_balance_subset, 
                                  chrX_pituitary$median_allele_balance_subset, 
                                  chrX_spleen$median_allele_balance_subset, 
                                  chrX_stomach$median_allele_balance_subset, 
                                  chrX_thyroid$median_allele_balance_subset, 
                                  chrX_uterus$median_allele_balance_subset, 
                                  chrX_vagina$median_allele_balance_subset, 
                                  chrX_placenta$median_allele_balance_subset), 
                       tissues = c(adrenal_lab,
                                   artery_lab,
                                   amygdala_lab,
                                   cortex_lab,
                                   breast_lab,
                                   colon_lab,
                                   heart_lab,
                                   kidney_lab,
                                   liver_lab,
                                   lung_lab,
                                   muscle_lab,
                                   nerve_lab,
                                   ovary_lab,
                                   pancreas_lab,
                                   pituitary_lab,
                                   spleen_lab,
                                   stomach_lab,
                                   thyroid_lab,
                                   uterus_lab,
                                   vagina_lab,
                                   placenta_lab),
                       type = c(gtex_lab, placenta_lab))

chrX_data$tissues = factor(chrX_data$tissues, levels = unique(chrX_data$tissues))

p1 = ggplot(chrX_data, aes(x=tissues, y=median, fill=type)) +
  geom_violin() +
  theme_bw() +
  coord_cartesian(ylim=c(0.5, 1)) +
  labs(y='Median allele balance', title = 'chrX') + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.title=element_blank(), legend.position="top", plot.title = element_text(hjust = 0.5, size=16), axis.title.x = element_blank(), legend.text = element_text(margin = margin(r = 30, unit = "pt"), size=16)) +
  geom_hline(yintercept = 0.5, linetype=2, col="black") +
  geom_hline(yintercept = 0.8, linetype=2, col="black")

# chr8
chr8_adrenal = read.csv('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chr8/AdrenalGland_chr8_median_allele_balance.txt', sep = '\t')
chr8_artery = read.csv('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chr8/ArteryAorta_chr8_median_allele_balance.txt', sep = '\t')
chr8_amygdala = read.csv('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chr8/BrainAmygdala_chr8_median_allele_balance.txt', sep = '\t')
chr8_cortex = read.csv('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chr8/BrainCortex_chr8_median_allele_balance.txt', sep = '\t')
chr8_breast = read.csv('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chr8/BreastMammaryTissue_chr8_median_allele_balance.txt', sep = '\t')
chr8_colon = read.csv('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chr8/ColonSigmoid_chr8_median_allele_balance.txt', sep = '\t')
chr8_heart = read.csv('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chr8/HeartLeftVentricle_chr8_median_allele_balance.txt', sep = '\t')
chr8_kidney = read.csv('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chr8/KidneyCortex_chr8_median_allele_balance.txt', sep = '\t')
chr8_liver = read.csv('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chr8/Liver_chr8_median_allele_balance.txt', sep = '\t')
chr8_lung = read.csv('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chr8/Lung_chr8_median_allele_balance.txt', sep = '\t')
chr8_muscle = read.csv('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chr8/MuscleSkeletal_chr8_median_allele_balance.txt', sep = '\t')
chr8_nerve = read.csv('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chr8/NerveTibial_chr8_median_allele_balance.txt', sep = '\t')
chr8_ovary = read.csv('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chr8/Ovary_chr8_median_allele_balance.txt', sep = '\t')
chr8_pancreas = read.csv('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chr8/Pancreas_chr8_median_allele_balance.txt', sep = '\t')
chr8_pituitary = read.csv('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chr8/Pituitary_chr8_median_allele_balance.txt', sep = '\t')
chr8_spleen = read.csv('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chr8/Spleen_chr8_median_allele_balance.txt', sep = '\t')
chr8_stomach = read.csv('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chr8/Stomach_chr8_median_allele_balance.txt', sep = '\t')
chr8_thyroid = read.csv('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chr8/Thyroid_chr8_median_allele_balance.txt', sep = '\t')
chr8_uterus = read.csv('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chr8/Uterus_chr8_median_allele_balance.txt', sep = '\t')
chr8_vagina = read.csv('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chr8/Vagina_chr8_median_allele_balance.txt', sep = '\t')
chr8_placenta = read.csv('/scratch/tphung3/PlacentaSexDiff/A_placenta/05_genotype_rnaseq/asereadcounter/HISAT/after_vqsr/chr8/placenta_chr8_median_allele_balance.txt', sep = '\t')

adrenal_lab = rep('Adrenal', 10)
artery_lab = rep('Artery', 10)
amygdala_lab = rep('Amygdala', 10)
cortex_lab = rep('Cortex', 10)
breast_lab = rep('Breast', 10)
colon_lab = rep('Colon', 10)
heart_lab = rep('Heart', 10)
kidney_lab = rep('Kidney', 10)
liver_lab = rep('Liver', 10)
lung_lab = rep('Lung', 10)
muscle_lab = rep('Muscle', 10)
nerve_lab = rep('Nerve', 10)
ovary_lab = rep('Ovary', 10)
pancreas_lab = rep('Pancreas', 10)
pituitary_lab = rep('Pituitary', 10)
spleen_lab = rep('Spleen', 10)
stomach_lab = rep('Stomach', 10)
thyroid_lab = rep('Thyroid', 10)
uterus_lab = rep('Uterus', 10)
vagina_lab = rep('Vagina', 10)
placenta_lab = rep('Placenta', 60)

gtex_lab = rep('GTEX', 200)

chr8_data = data.frame(median = c(chr8_adrenal$median_allele_balance_subset, 
                                  chr8_artery$median_allele_balance_subset, 
                                  chr8_amygdala$median_allele_balance_subset, 
                                  chr8_cortex$median_allele_balance_subset, 
                                  chr8_breast$median_allele_balance_subset, 
                                  chr8_colon$median_allele_balance_subset, 
                                  chr8_heart$median_allele_balance_subset, 
                                  chr8_kidney$median_allele_balance_subset, 
                                  chr8_liver$median_allele_balance_subset, 
                                  chr8_lung$median_allele_balance_subset, 
                                  chr8_muscle$median_allele_balance_subset, 
                                  chr8_nerve$median_allele_balance_subset, 
                                  chr8_ovary$median_allele_balance_subset, 
                                  chr8_pancreas$median_allele_balance_subset, 
                                  chr8_pituitary$median_allele_balance_subset, 
                                  chr8_spleen$median_allele_balance_subset, 
                                  chr8_stomach$median_allele_balance_subset, 
                                  chr8_thyroid$median_allele_balance_subset, 
                                  chr8_uterus$median_allele_balance_subset, 
                                  chr8_vagina$median_allele_balance_subset, 
                                  chr8_placenta$median_allele_balance_subset), 
                       tissues = c(adrenal_lab,
                                   artery_lab,
                                   amygdala_lab,
                                   cortex_lab,
                                   breast_lab,
                                   colon_lab,
                                   heart_lab,
                                   kidney_lab,
                                   liver_lab,
                                   lung_lab,
                                   muscle_lab,
                                   nerve_lab,
                                   ovary_lab,
                                   pancreas_lab,
                                   pituitary_lab,
                                   spleen_lab,
                                   stomach_lab,
                                   thyroid_lab,
                                   uterus_lab,
                                   vagina_lab,
                                   placenta_lab),
                       type = c(gtex_lab, placenta_lab))

chr8_data$tissues = factor(chr8_data$tissues, levels = unique(chr8_data$tissues))

p2 = ggplot(chr8_data, aes(x=tissues, y=median, fill=type)) +
  geom_violin() +
  theme_bw() +
  coord_cartesian(ylim=c(0.5, 1)) +
  labs(y='Median allele balance', title = 'chr8') + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.title=element_blank(), legend.position="top", plot.title = element_text(hjust = 0.5, size=16), axis.title.x = element_blank(), legend.text = element_text(margin = margin(r = 30, unit = "pt"), size=16)) +
  geom_hline(yintercept = 0.5, linetype=2, col="black") +
  geom_hline(yintercept = 0.8, linetype=2, col="black")

png('/scratch/tphung3/PlacentaSexDiff/B_gtex/08_genotype_rnaseq/compare_violin_gtex_placenta_rnaseq_genotypes.png', width = 15, height = 9, units = "in", res = 300)
ggarrange(p2, p1, ncol = 1, nrow = 2, common.legend = T)
dev.off()