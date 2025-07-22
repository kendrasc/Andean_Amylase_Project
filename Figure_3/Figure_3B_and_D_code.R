library(ggplot2)
library(dplyr)
library(ggpubr)


setwd("~/Desktop/Amylase_Americas/")
Genotype = read.delim("AllQuechua_Maya_PhasedChr1_no_missing_biallelic_no_EUR_maf_0.05_genotypes_recreate.txt", header = T)
GeneToCN = read.delim("All_amylase_gene_copy_number.txt", header = T)

# filter to just include SNPs in the amylase region
Genotype = Genotype %>%
  dplyr::filter(POS > 103343000 & POS < 103960000 )

# organize genotype data
Genotype = subset(Genotype, select = -CHROM)
Genotype = Genotype %>%
  mutate(POS = paste0("X", POS))
rownames(Genotype) = Genotype$POS
Genotype_transposed = as.data.frame(t(Genotype[,-1])) 
colnames(Genotype_transposed) = Genotype$POS
Genotype_transposed$ID = rownames(Genotype_transposed) 
rownames(Genotype_transposed) = NULL
Genotype_transposed = Genotype_transposed %>%
  select(ID, everything())

# Organize populations
GenetoCN_Pops = GeneToCN %>%
  dplyr::filter(Pop %in% c("Lowland_Andeans", "Highland_Andeans", "Maya_Chiapas")) %>%
  dplyr::mutate(Pop = dplyr::case_when(
    Pop %in% c("Lowland_Andeans", "Highland_Andeans") ~ "Quechua",
    TRUE ~ Pop
  ))
################################# SNP by Gene Copy Number ###############################
GenetoCN_Pops_sorted = GenetoCN_Pops[order(GenetoCN_Pops$ID),]

Genotype_transposed = Genotype_transposed %>% 
  mutate(ID_short = sub(".*_", "", ID)) %>% 
  select(ID_short, everything())

Genotype_sorted = Genotype_transposed[order(Genotype_transposed$ID_short),]
dataset_for_snp_analyses = cbind(GenetoCN_Pops_sorted, Genotype_sorted[, -2])

dataset_for_snp_analyses[9:ncol(dataset_for_snp_analyses)] = lapply(
  dataset_for_snp_analyses[9:ncol(dataset_for_snp_analyses)],
  function(x) factor(ifelse(x == "0|0", "Homozygous Ancestral",
                            ifelse(x == "1|1", "Homozygous Derived", "Heterozygous")),
                     levels = c("Homozygous Ancestral", "Heterozygous", "Homozygous Derived")))


###################### Top SNP by AMY1 Copy Number ##################
my_comparisons = list( c("Homozygous Ancestral", "Heterozygous"), c("Heterozygous", "Homozygous Derived"), c("Homozygous Ancestral", "Homozygous Derived") )

  # Example SNV is the SNV from Figures 3B and 3C
figure = dataset_for_snp_analyses %>%
  ggplot(., aes(x = X103614521, y = round(AMY1), color = Pop)) +
  geom_boxplot() +
  stat_compare_means(comparisons = my_comparisons) +
  scale_color_manual(values = rep(c("#800080", "#D4AA00"), 22 )) +
  labs(color = "Population") +
  ylim(0, 29) +
  ggtitle("AMY1 Copy number by SNP") +
  xlab("SNV Position 103614521") +
  ylab("AMY1 Copy Number") +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75)) +
  theme_classic() +
  theme( legend.position = "bottom")

figure
pdf(file = "SNP_figure.pdf", width = 5, height = 4)
print(figure)
dev.off()

