library(dplyr)
library(data.table)
library(ggplot2)
library(viridis)

FST = read.delim("Quechua_Maya_amylase_all_chr1_biallelic_0.05_1_region_fst.txt", header = T)
Genotype = read.delim("1k_genomes_Quechua_Maya_chr1_no_indels_0.05_no_missing_unique_Genotype.txt", header = T)
GeneToCN = read.delim("All_amylase_gene_copy_number.txt", header = T)

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

saveRDS(Genotype_transposed, "AMR_chr1_snps_Quechua_Maya_reordered.rds")

GenetoCN_Pops = GeneToCN %>%
  dplyr::filter(Pop == "Maya" | Pop == "Quechua")

################################# SNP by Gene Copy Number ###############################
GenetoCN_Pops_sorted = GenetoCN_Pops[order(GenetoCN_Pops$ID),]
Genotype_sorted = Genotype_transposed[order(Genotype_transposed$ID),]
dataset_for_snp_analyses = cbind(GenetoCN_Pops_sorted, Genotype_sorted[, 2:ncol(Genotype_sorted)])
snp_data <- dataset_for_snp_analyses[, 8:length(dataset_for_snp_analyses)]

saveRDS(dataset_for_snp_analyses, "AMR_chr1_snps_Quechua_Maya_reordered_with_copy_number.rds")

################### Relabel SNPs to combined Heterozygous Individuals ###############################
dataset_for_snp_analyses[8:length(dataset_for_snp_analyses)] <- lapply(
  snp_data,
  function(x) factor(
    ifelse(x == "0|0", "Homozygous Ancestral",
           ifelse(x == "1|1", "Homozygous Derived", "Heterozygous")),
    levels = c("Homozygous Ancestral", "Heterozygous", "Homozygous Derived")
  )
)

################### Save in case of later issues ####################
saveRDS(dataset_for_snp_analyses, "AMR_chr1_snps_Quechua_Maya_reordered_with_copy_number_relabeled.rds")

################### Calculating Kruskal Wallis p-value between gene copy number and snp ###############
AMY1_data <- dataset_for_snp_analyses$AMY1

Kruskal_pvalue <- vapply(snp_data, function(snp) {
  kruskal.test(AMY1_data ~ snp)$p.value
}, numeric(1))

############ Check to make sure numeric ##############
Kruskal_pvalue <- as.numeric(Kruskal_pvalue)

##################### X axis #######################
Kruskal_Wallis_status = numeric()
Kruskal_Wallis_status <- sapply(dataset_for_snp_analyses[, 8:length(dataset_for_snp_analyses)], function(col) {
  
  # Filter values in gene copy number column based on current column
  Hets <- dataset_for_snp_analyses[[3]][col == "Heterozygous"]
  Ancestral <- dataset_for_snp_analyses[[3]][col == "Homozygous Ancestral"]
  Derived <- dataset_for_snp_analyses[[3]][col == "Homozygous Derived"]
  # Calculate medians for the filtered groups, handling empty groups
  Het_median <- if (length(Hets) > 0) median(Hets) else NA
  Ancestral_median <- if (length(Ancestral) > 0) median(Ancestral) else NA
  Derived_median <- if (length(Derived) > 0) median(Derived) else NA
  
  # Describe conditions based on medians, make sure there aren't any NAs
  if (!is.na(Derived_median) && !is.na(Het_median) && !is.na(Ancestral_median)) {
    condition_pos_met <- Derived_median > Het_median && Het_median > Ancestral_median
    condition_neg_met <- Derived_median < Het_median && Het_median < Ancestral_median
  } else if (!is.na(Derived_median) && !is.na(Het_median) && is.na(Ancestral_median)) {
    condition_pos_met <- Derived_median > Het_median
    condition_neg_met <- Derived_median < Het_median
  } else if (is.na(Derived_median) && !is.na(Het_median) && !is.na(Ancestral_median)) {
    condition_pos_met <- Het_median > Ancestral_median
    condition_neg_met <- Het_median < Ancestral_median
  } else if (!is.na(Derived_median) && is.na(Het_median) && !is.na(Ancestral_median)) {
    condition_pos_met <- Derived_median > Ancestral_median
    condition_neg_met <- Derived_median < Ancestral_median
  }
  else {
    # If any median is NA, no valid comparison can be made
    return(0)
  }
  
  # Return value based on conditions
  if (condition_pos_met) {
    return(1)
  } else if (condition_neg_met) {
    return(-1)
  } else {
    return(0)
  }
})

##################### Y axis ########################
Quechua <- as.data.table(Quechua)
Maya <- as.data.table(Maya)

# Determine FST direct (positive to first population, here Quechua, negative to second population, here Maya)
FST_status <- sapply(names(dataset_for_snp_analyses)[8:length(dataset_for_snp_analyses)], function(col_name) {
  
  # Extract columns once
  col_Quechua <- Quechua[[col_name]]
  col_Maya <- Maya[[col_name]]
  
  # Use table() to count all values at once
  count_Quechua <- table(factor(col_Quechua, levels = c("Heterozygous", "Homozygous Ancestral", "Homozygous Derived")))
  count_Maya <- table(factor(col_Maya, levels = c("Heterozygous", "Homozygous Ancestral", "Homozygous Derived")))
  
  # Extract counts
  Hets_Quechua <- count_Quechua["Heterozygous"]
  Ancestral_Quechua <- count_Quechua["Homozygous Ancestral"]
  Derived_Quechua <- count_Quechua["Homozygous Derived"]
  
  Hets_Maya <- count_Maya["Heterozygous"]
  Ancestral_Maya <- count_Maya["Homozygous Ancestral"]
  Derived_Maya <- count_Maya["Homozygous Derived"]
  
  # Replace NAs with 0 for count values
  Hets_Quechua <- ifelse(is.na(Hets_Quechua), 0, Hets_Quechua)
  Ancestral_Quechua <- ifelse(is.na(Ancestral_Quechua), 0, Ancestral_Quechua)
  Derived_Quechua <- ifelse(is.na(Derived_Quechua), 0, Derived_Quechua)
  
  Hets_Maya <- ifelse(is.na(Hets_Maya), 0, Hets_Maya)
  Ancestral_Maya <- ifelse(is.na(Ancestral_Maya), 0, Ancestral_Maya)
  Derived_Maya <- ifelse(is.na(Derived_Maya), 0, Derived_Maya)
  
  # Calculate ancestral frequencies
  total_Quechua <- Hets_Quechua + Ancestral_Quechua + Derived_Quechua
  total_Maya <- Hets_Maya + Ancestral_Maya + Derived_Maya
  
  if (total_Quechua == 0 || total_Maya == 0) {
    return(0)
  }
  
  Frequency_Quechua_ancestral <- (Hets_Quechua + 2 * Derived_Quechua) / (2 * total_Quechua)
  Frequency_Maya_ancestral <- (Hets_Maya + 2 * Derived_Maya) / (2 * total_Maya)
  
  # Determine FST status
  if (Frequency_Quechua_ancestral > Frequency_Maya_ancestral) {
    1
  } else if (Frequency_Maya_ancestral > Frequency_Quechua_ancestral) {
    -1
  } else {
    0
  }
})

###################### P-value by FST graph #######################
FST$WEIR_AND_COCKERHAM_FST[FST$WEIR_AND_COCKERHAM_FST < 0] <- 0
FST_directional = FST$WEIR_AND_COCKERHAM_FST*FST_status
Kruskal_Wallis_status = as.numeric(unlist(Kruskal_Wallis_status))
Krusal_pvalue_directional = (-log10(Kruskal_pvalue)*Kruskal_Wallis_status)
FST_by_KW_test = cbind(FST[1:2], Krusal_pvalue_directional, FST_directional)
FST_by_KW_test$Color_the_figure = FST_by_KW_test$Krusal_pvalue_directional * FST_by_KW_test$FST_directional


figure = FST_by_KW_test %>%
  dplyr::filter(Krusal_pvalue_directional != 0 & FST_directional != 0) %>%
  ggplot(., aes(x = Krusal_pvalue_directional, y = FST_directional, color = Color_the_figure)) +
  geom_point() + 
  theme_classic() +
  scale_color_viridis_c(direction = -1, option = "inferno") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0.1265561, linetype = "dotted") +
  geom_hline(yintercept = -0.1265561, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  xlab("-log(Kruskal-Wallis p-value) SNP by AMY1 Copy") +
  ylab("SNP FST") +
  xlim(-18, 18) +
  ylim(-0.37, 0.37) +
  labs(color = "FST by P-value") +
  ggtitle("Quechua vs. Maya Derived SNP Status: All Chr. 1 SNPs")
figure

pdf(file = "Figure_3_Quechua_Maya.pdf", width = 10, height = 8)
print(figure)
dev.off()


