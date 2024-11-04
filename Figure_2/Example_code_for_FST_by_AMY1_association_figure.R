# Code to Create Figure for Supplementary Figure X  Figure 2A

FST = read.delim("KHV_JPT_amylase_all_chr1_biallelic_0.05_1_region_fst.txt", header = T)
Genotype = read.delim("1k_genomes_KHV_JPT_chr1_no_indels_0.05_no_missing_unique_Genotype.txt", header = T)
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

saveRDS(Genotype_transposed, "EAS_chr1_snps_KHV_JPT_reordered.rds")

GenetoCN_Pops = GeneToCN %>%
  dplyr::filter(Pop == "JPT" | Pop == "KHV")

################################# SNP by Gene Copy Number ###############################
GenetoCN_Pops_sorted = GenetoCN_Pops[order(GenetoCN_Pops$ID),]
Genotype_sorted = Genotype_transposed[order(Genotype_transposed$ID),]
dataset_for_snp_analyses = cbind(GenetoCN_Pops_sorted, Genotype_sorted[, 2:ncol(Genotype_sorted)])
snp_data <- dataset_for_snp_analyses[, 8:length(dataset_for_snp_analyses)]

saveRDS(dataset_for_snp_analyses, "EAS_chr1_snps_KHV_JPT_reordered_with_copy_number.rds")

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
saveRDS(dataset_for_snp_analyses, "EAS_chr1_snps_KHV_JPT_reordered_with_copy_number_relabeled.rds")

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
library(data.table)

KHV <- as.data.table(KHV)
JPT <- as.data.table(JPT)

# Determine FST direct (positive to first population, here KHV, negative to second population, here JPT)
FST_status <- sapply(names(dataset_for_snp_analyses)[8:length(dataset_for_snp_analyses)], function(col_name) {
  
  # Extract columns once
  col_KHV <- KHV[[col_name]]
  col_JPT <- JPT[[col_name]]
  
  # Use table() to count all values at once
  count_KHV <- table(factor(col_KHV, levels = c("Heterozygous", "Homozygous Ancestral", "Homozygous Derived")))
  count_JPT <- table(factor(col_JPT, levels = c("Heterozygous", "Homozygous Ancestral", "Homozygous Derived")))
  
  # Extract counts
  Hets_KHV <- count_KHV["Heterozygous"]
  Ancestral_KHV <- count_KHV["Homozygous Ancestral"]
  Derived_KHV <- count_KHV["Homozygous Derived"]
  
  Hets_JPT <- count_JPT["Heterozygous"]
  Ancestral_JPT <- count_JPT["Homozygous Ancestral"]
  Derived_JPT <- count_JPT["Homozygous Derived"]
  
  # Replace NAs with 0 for count values
  Hets_KHV <- ifelse(is.na(Hets_KHV), 0, Hets_KHV)
  Ancestral_KHV <- ifelse(is.na(Ancestral_KHV), 0, Ancestral_KHV)
  Derived_KHV <- ifelse(is.na(Derived_KHV), 0, Derived_KHV)
  
  Hets_JPT <- ifelse(is.na(Hets_JPT), 0, Hets_JPT)
  Ancestral_JPT <- ifelse(is.na(Ancestral_JPT), 0, Ancestral_JPT)
  Derived_JPT <- ifelse(is.na(Derived_JPT), 0, Derived_JPT)
  
  # Calculate ancestral frequencies
  total_KHV <- Hets_KHV + Ancestral_KHV + Derived_KHV
  total_JPT <- Hets_JPT + Ancestral_JPT + Derived_JPT
  
  if (total_KHV == 0 || total_JPT == 0) {
    return(0) # Prevent division by zero if no data is available
  }
  
  Frequency_KHV_ancestral <- (Hets_KHV + 2 * Derived_KHV) / (2 * total_KHV)
  Frequency_JPT_ancestral <- (Hets_JPT + 2 * Derived_JPT) / (2 * total_JPT)
  
  # Determine FST status
  if (Frequency_KHV_ancestral > Frequency_JPT_ancestral) {
    1
  } else if (Frequency_JPT_ancestral > Frequency_KHV_ancestral) {
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

################ 
saveRDS(FST_by_KW_test, "EAS_chr1_snps_KHV_JPT_reordered_with_copy_number_relabeled_Final_dataset.rds")

figure = FST_by_KW_test %>%
  dplyr::filter(Krusal_pvalue_directional != 0 & FST_directional != 0) %>%
  ggplot(., aes(x = Krusal_pvalue_directional, y = FST_directional, color = Color_the_figure)) +
  geom_point() + 
  theme_classic() +
  scale_color_viridis_c(direction = -1, option = "inferno") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  xlab("-log(Kruskal-Wallis p-value) SNP by AMY1 Copy") +
  ylab("SNP FST") +
  xlim(-18, 18) +
  ylim(-0.42, 0.42) +
  labs(color = "FST by P-value") +
  ggtitle("KHV vs. JPT Derived SNP Status: All Chr. 1 SNPs")

pdf(file = "Cross_figure_KHV_JPT.pdf", width = 10, height = 8)
print(figure)
dev.off()


