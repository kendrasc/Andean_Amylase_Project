library(ggplot2)
library(dplyr)
library(ggpubr)

setwd("~/Desktop/Amylase_Americas/")
Genotype = read.delim("AllQuechua_Maya_PhasedChr1_no_missing_biallelic_no_EUR_maf_0.05_genotypes.txt", header = T)
GeneToCN = read.delim("All_amylase_gene_copy_number.txt", header = T)

Genotype = Genotype %>%
  dplyr::filter(POS > 103343000 & POS < 103960000 )

Genotype <- subset(Genotype, select = -CHROM)
Genotype <- Genotype %>%
  mutate(POS = paste0("X", POS))
rownames(Genotype) = Genotype$POS
Genotype_transposed <- as.data.frame(t(Genotype[,-1])) 
colnames(Genotype_transposed) <- Genotype$POS
Genotype_transposed$ID <- rownames(Genotype_transposed) 
rownames(Genotype_transposed) <- NULL
Genotype_transposed <- Genotype_transposed %>%
  select(ID, everything())

GenetoCN_Pops <- GeneToCN %>%
  dplyr::filter(Pop %in% c("Lowland_Andeans", "Highland_Andeans", "Maya_Chiapas")) %>%
  # Remove Quechua individuals with 100% non-AMR ancestry here
  dplyr::mutate(Pop = dplyr::case_when(
    Pop %in% c("Lowland_Andeans", "Highland_Andeans") ~ "Quechua",
    TRUE ~ Pop
  ))
################################# SNP by Gene Copy Number ###############################
GenetoCN_Pops_sorted = GenetoCN_Pops[order(GenetoCN_Pops$ID),]

Genotype_transposed <- Genotype_transposed %>% 
  mutate(ID_short = sub(".*_", "", ID)) %>% 
  select(ID_short, everything())

Genotype_sorted = Genotype_transposed[order(Genotype_transposed$ID_short),]
dataset_for_snp_analyses <- cbind(GenetoCN_Pops_sorted, Genotype_sorted[, -2])

dataset_for_snp_analyses[9:ncol(dataset_for_snp_analyses)] <- lapply(
  dataset_for_snp_analyses[9:ncol(dataset_for_snp_analyses)],
  function(x) factor(ifelse(x == "0|0", "Homozygous Ancestral",
                            ifelse(x == "1|1", "Homozygous Derived", "Heterozygous")),
                     levels = c("Homozygous Ancestral", "Heterozygous", "Homozygous Derived")))


############################## Correlation testing
results <- lapply(names(dataset_for_snp_analyses)[9:ncol(dataset_for_snp_analyses)], function(col_name) {
  x <- as.numeric(dataset_for_snp_analyses[[col_name]])
  y <- round(dataset_for_snp_analyses[[3]])
  
  res <- cor.test(x, y, method = "pearson")
  data.frame(SNP = col_name, Correlation = res$estimate, P = res$p.value)
})

cor_df <- do.call(rbind, results)
cor_df_sorted <- cor_df[order(cor_df$Correlation), ]


################ Quechua as an example
Quechua = dataset_for_snp_analyses %>%
  dplyr::filter(Pop == "Quechua")

m_results <- lapply(names(Quechua)[9:ncol(Quechua)], function(col_name) {
  x <- as.numeric(Quechua[[col_name]])
  y <- round(Quechua[[3]])
  
  res <- cor.test(x, y, method = "pearson")
  data.frame(SNP = col_name, Correlation = res$estimate, P = res$p.value)
})

m_cor_df <- do.call(rbind, m_results)
m_cor_df_sorted <- m_cor_df[order(m_cor_df$Correlation), ]

#example SNV
m_cor_df_sorted %>%
  dplyr::filter(SNP == "X103614521")
