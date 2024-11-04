library(ggplot2)
library(dplyr)
library(ggpubr)

# Removed the last position from FST b/c part of another haplotype and for some reason, while it says it is being used in the log file, it doesn't appear in the LD analysis
Genotype = read.delim("Quechua_Missing_EUR_in_amylase_Maya_103343_10396_biallelic_0.05_1_region_Genotypes.txt", header = T)
GeneToCN = read.delim("All_Amylase_I_think.txt", header = T)

GenetoCN_Indigenous = GeneToCN %>%
  dplyr::filter(Pop == "Lowland_Andeans" | Pop == "Highland_Andeans" | Pop == "Maya_Chiapas") %>%
  dplyr::filter (ID != "108333" & ID != "108337" & ID != "CNS0261883")

Andean = c("Lowland_Andeans", "Highland_Andeans")
GenetoCN_Indigenous$Andean = with(GenetoCN_Indigenous,
                                       ifelse(Pop == "Lowland_Andeans" | Pop == "Highland_Andeans", "Andean", "Maya"))
Andeans = as.data.frame(cbind(GenetoCN_Indigenous$ID, GenetoCN_Indigenous$AMY1, GenetoCN_Indigenous$Andean))

################################# SNP by Gene Copy Number ###############################
GenetoCN_Indigenous_sorted = GenetoCN_Indigenous[order(GenetoCN_Indigenous$ID),]
Genotype_sorted = Genotype[order(Genotype$ID),]
dataset_for_snp_analyses = cbind(GenetoCN_Indigenous_sorted, Genotype_sorted[, 2:ncol(Genotype_sorted)])

dataset_for_snp_analyses[8:length(dataset_for_snp_analyses)] <- lapply(
  dataset_for_snp_analyses[8:length(dataset_for_snp_analyses)],
  function(x) factor(
    ifelse(x == "0|0", "Homozygous Ancestral",
           ifelse(x == "1|1", "Homozygous Derived", "Heterozygous")),
    levels = c("Homozygous Ancestral", "Heterozygous", "Homozygous Derived")
  )
)

################### Calculating Kruskal Wallis p-value between gene copy number and snp ###############
Andean = c("Lowland_Andeans", "Highland_Andeans")
dataset_for_snp_analyses$Andean = with(dataset_for_snp_analyses,
                                       ifelse(Pop == "Lowland_Andeans" | Pop == "Highland_Andeans", "Andean", "Maya"))
###################### Top SNP by AMY1 Copy Number ##################
my_comparisons <- list( c("Homozygous Ancestral", "Heterozygous"), c("Heterozygous", "Homozygous Derived"), c("Homozygous Ancestral", "Homozygous Derived") )
my_comparisons_2 <- c("Andean", "Maya")


figure = dataset_for_snp_analyses %>%
  ggplot(., aes(x = X103584693  , y = round(AMY1), color = Andean)) +
  geom_boxplot() +
  stat_compare_means(comparisons = my_comparisons) +
  stat_compare_means(aes(group = Andean), method = "wilcox.test") +
  scale_color_manual(values = rep(c("orange", "maroon"), 22 )) +
  labs(color = "Population") +
  ylim(0, 29) +
  ggtitle("AMY1 Copy number by SNP") +
  xlab("SNP Position 103584693") +
  ylab("AMY1 Copy Number") +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75)) +
  theme_classic()


pdf(file = "103584693_SNP_figure.pdf", width = 6, height = 4)
print(figure)
dev.off()






