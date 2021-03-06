library(ggplot2)
library(ggpubr)

args = commandArgs(trailingOnly=TRUE)
tissue = args[1]
# tissue = 'AdrenalGland'
# ----
# chr8
# ----
chr8 = read.csv(paste('/scratch/tphung3/PlacentaSexDiff/B_gtex/09_compare_wes_vs_rnaseq_genotypes/results/chr8/compare_wes_rna_genotyping_chr8_', tissue, '.txt', sep = ''), sep='\t')

head(chr8)

chr8_df = data.frame(median = c(
  chr8$wes_subset_median,
  chr8$rna_subset_median,
  chr8$wes_subset_shared_median,
  chr8$rna_subset_shared_median,
  chr8$wes_subset_unique_median,
  chr8$rna_subset_unique_median
),
categories = c(rep('All', 20), rep('Shared', 20), rep('Unique', 20)),
labels = c(rep('WES_All', 10), rep('RNA_All', 10), rep('WES_Shared', 10), rep('RNA_Shared', 10), rep('WES_Unique', 10), rep('RNA_Unique', 10))
)

chr8_df$labels = factor(chr8_df$labels, levels = unique(chr8_df$labels))

p1 = ggplot(chr8_df, aes(x=labels, y=median, fill=categories)) +
  geom_violin() +
  theme_bw() +
  coord_cartesian(ylim=c(0.5, 1)) +
  labs(y='Median allele balance', title = paste('chr8_', tissue, sep='')) + 
  theme(legend.title=element_blank(), legend.position="top", plot.title = element_text(hjust = 0.5, size=20), axis.title.x = element_blank(), legend.text = element_text(margin = margin(r = 30, unit = "pt"), size=16), axis.text = element_text(size=16)) +
  geom_hline(yintercept = 0.5, linetype=2, col="black") +
  geom_hline(yintercept = 0.8, linetype=2, col="black")

# ----
# chrX
# ----
chrX = read.csv(paste('/scratch/tphung3/PlacentaSexDiff/B_gtex/09_compare_wes_vs_rnaseq_genotypes/results/chrX/compare_wes_rna_genotyping_chrX_', tissue, '.txt', sep = ''), sep='\t')

# head(chrX)

chrX_df = data.frame(median = c(
  chrX$wes_subset_median,
  chrX$rna_subset_median,
  chrX$wes_subset_shared_median,
  chrX$rna_subset_shared_median,
  chrX$wes_subset_unique_median,
  chrX$rna_subset_unique_median
),
categories = c(rep('All', 20), rep('Shared', 20), rep('Unique', 20)),
labels = c(rep('WES_All', 10), rep('RNA_All', 10), rep('WES_Shared', 10), rep('RNA_Shared', 10), rep('WES_Unique', 10), rep('RNA_Unique', 10))
)

chrX_df$labels = factor(chrX_df$labels, levels = unique(chrX_df$labels))

p2 = ggplot(chrX_df, aes(x=labels, y=median, fill=categories)) +
  geom_violin() +
  theme_bw() +
  coord_cartesian(ylim=c(0.5, 1)) +
  labs(y='Median allele balance', title = paste('chrX_', tissue, sep='')) + 
  theme(legend.title=element_blank(), legend.position="top", plot.title = element_text(hjust = 0.5, size=20), axis.title.x = element_blank(), legend.text = element_text(margin = margin(r = 30, unit = "pt"), size=16), axis.text = element_text(size=16)) +
  geom_hline(yintercept = 0.5, linetype=2, col="black") +
  geom_hline(yintercept = 0.8, linetype=2, col="black")

png(paste('/scratch/tphung3/PlacentaSexDiff/B_gtex/09_compare_wes_vs_rnaseq_genotypes/results/compare_genotype_calling_from_wes_rnaseq_', tissue, '.png', sep=''), width = 15, height = 9, units = "in", res = 300)
ggarrange(p1, p2, ncol = 1, nrow = 2, common.legend = T)
dev.off()