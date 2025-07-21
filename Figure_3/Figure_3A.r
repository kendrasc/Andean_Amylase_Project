library(ggplot2)
library(dplyr)
library(ggpubr)
library(forcats)
library(FSA)
library(data.table)

setwd("~/Desktop/Amylase_Americas/")

# There are many saves within this code since the script kept crashing/having issues while I was writing it
# This code could also probably be a lot cleaner but I'm mostly just meh at coding...

FST = read.delim("Quechua_Missing_EUR_in_amylase_Maya_all_chr1_biallelic_0.05_1_region_fst.txt", header = T)
Genotype = read.delim("AllQuechua_Maya_PhasedChr1_no_missing_biallelic_no_EUR_maf_0.05_genotypes_recreated.txt", header = T)
GeneToCN = read.delim("All_amylase_gene_copy_number.txt", header = T)

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

saveRDS(Genotype_transposed, "AMR_chr1_snps_Quechua_Maya_reordered.rds")

GenetoCN_Pops = GeneToCN %>%
  dplyr::filter(Pop %in% c("Lowland_Andeans", "Highland_Andeans", "Maya_Chiapas")) %>%
  dplyr::filter(ID != "CNS0261883" & ID != "108337" & ID != "108333") %>%
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
dataset_for_snp_analyses = cbind(GenetoCN_Pops_sorted, Genotype_sorted[, -1])

dataset_for_snp_analyses[9:ncol(dataset_for_snp_analyses)] = lapply(
  dataset_for_snp_analyses[9:ncol(dataset_for_snp_analyses)],
  function(x) factor(ifelse(x == "0|0", "Homozygous Ancestral",
                            ifelse(x == "1|1", "Homozygous Derived", "Heterozygous")),
                     levels = c("Homozygous Ancestral", "Heterozygous", "Homozygous Derived"))
)

saveRDS(dataset_for_snp_analyses, "AMR_chr1_snps_Quechua_Maya_reordered_with_copy_number.rds")
dataset_for_snp_analyses$AMY1 = round(dataset_for_snp_analyses$AMY1)

Quechua_indices = dataset_for_snp_analyses[[2]] == "Quechua"
Maya_indices = dataset_for_snp_analyses[[2]] == "Maya_Chiapas"
Quechua = dataset_for_snp_analyses[Quechua_indices, ]
Maya = dataset_for_snp_analyses[Maya_indices, ]

AMY1_data = dataset_for_snp_analyses$AMY1
snp_data = dataset_for_snp_analyses[, 9:length(dataset_for_snp_analyses)]

Kruskal_pvalue = vapply(snp_data, function(snp) {
  kruskal.test(AMY1_data ~ snp)$p.value
}, numeric(1))

saveRDS(Kruskal_pvalue, "AMR_chr1_snps_Quechua_Maya_Kruskal_pvalue.rds")

KW_adj = p.adjust(Kruskal_pvalue, method = "bonferroni") 
write.csv(KW_adj, "Kruskal_Wallis_Bonferroni_adjusted_p_values.csv")

sig_snps = names(KW_adj)[KW_adj < 0.05]

dunn_list = lapply(sig_snps, function(snp) {
  res <- dunnTest(AMY1_data ~ dataset_for_snp_analyses[[snp]],
                  method = "bonferroni")$res
  cbind(SNP = snp,
        KW_q = KW_adj[[snp]],
        res)
})

dunn_results = do.call(rbind, dunn_list)
rownames(dunn_results) = NULL   
write.csv(dunn_results, "Dunn_test_Quechua_Maya_results.csv")

Kruskal_pvalue = as.numeric(Kruskal_pvalue)

##################### X axis #######################
Kruskal_Wallis_status = numeric()
Kruskal_Wallis_status = sapply(dataset_for_snp_analyses[, 9:length(dataset_for_snp_analyses)], function(col) {
  
  Hets = dataset_for_snp_analyses[[3]][col == "Heterozygous"]
  Ancestral = dataset_for_snp_analyses[[3]][col == "Homozygous Ancestral"]
  Derived = dataset_for_snp_analyses[[3]][col == "Homozygous Derived"]
  
  Het_median = if (length(Hets) > 0) median(Hets) else NA
  Ancestral_median = if (length(Ancestral) > 0) median(Ancestral) else NA
  Derived_median = if (length(Derived) > 0) median(Derived) else NA
  
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
    return(0)
  }
  
  if (condition_pos_met) {
    return(1)
  } else if (condition_neg_met) {
    return(-1)
  } else {
    return(0)
  }
})

saveRDS(Kruskal_Wallis_status, "AMR_chr1_snps_Quechua_Maya_Kruskal_Wallis_status.rds")

##################### Y axis testing #######################
Quechua = as.data.table(Quechua)
Maya = as.data.table(Maya)

FST_status = sapply(names(dataset_for_snp_analyses)[9:length(dataset_for_snp_analyses)], function(col_name) {
  
  col_Quechua = Quechua[[col_name]]
  col_Maya = Maya[[col_name]]
  
  count_Quechua = table(factor(col_Quechua, levels = c("Heterozygous", "Homozygous Ancestral", "Homozygous Derived")))
  count_Maya = table(factor(col_Maya, levels = c("Heterozygous", "Homozygous Ancestral", "Homozygous Derived")))
  
  Hets_Quechua = count_Quechua["Heterozygous"]
  Ancestral_Quechua = count_Quechua["Homozygous Ancestral"]
  Derived_Quechua = count_Quechua["Homozygous Derived"]
  
  Hets_Maya = count_Maya["Heterozygous"]
  Ancestral_Maya = count_Maya["Homozygous Ancestral"]
  Derived_Maya = count_Maya["Homozygous Derived"]
  
  Hets_Quechua = ifelse(is.na(Hets_Quechua), 0, Hets_Quechua)
  Ancestral_Quechua = ifelse(is.na(Ancestral_Quechua), 0, Ancestral_Quechua)
  Derived_Quechua = ifelse(is.na(Derived_Quechua), 0, Derived_Quechua)
  
  Hets_Maya = ifelse(is.na(Hets_Maya), 0, Hets_Maya)
  Ancestral_Maya = ifelse(is.na(Ancestral_Maya), 0, Ancestral_Maya)
  Derived_Maya = ifelse(is.na(Derived_Maya), 0, Derived_Maya)
  
  total_Quechua = Hets_Quechua + Ancestral_Quechua + Derived_Quechua
  total_Maya = Hets_Maya + Ancestral_Maya + Derived_Maya
  
  if (total_Quechua == 0 || total_Maya == 0) {
    return(0)}
  
  Frequency_Quechua_ancestral = (Hets_Quechua + 2 * Derived_Quechua) / (2 * total_Quechua)
  Frequency_Maya_ancestral = (Hets_Maya + 2 * Derived_Maya) / (2 * total_Maya)
  
  if (Frequency_Quechua_ancestral > Frequency_Maya_ancestral) {
    1
  } else if (Frequency_Maya_ancestral > Frequency_Quechua_ancestral) {
    -1
  } else {
    0
  }
})


saveRDS(FST_status, "AMR_chr1_snps_Quechua_Maya_FST_status.rds")

###################### P-value by FST graph #######################
FST$WEIR_AND_COCKERHAM_FST[FST$WEIR_AND_COCKERHAM_FST < 0] = 0
FST_directional = FST$WEIR_AND_COCKERHAM_FST*FST_status
Kruskal_Wallis_status = as.numeric(unlist(Kruskal_Wallis_status))
Krusal_pvalue_directional = (-log10(Kruskal_pvalue)*Kruskal_Wallis_status)
FST_by_KW_test = cbind(FST[1:2], Krusal_pvalue_directional, FST_directional)
FST_by_KW_test$Color_the_figure = FST_by_KW_test$Krusal_pvalue_directional * FST_by_KW_test$FST_directional
saveRDS(FST_by_KW_test, "AMR_chr1_snps_Quechua_Maya_reordered_with_copy_number_relabeled_Final_dataset.rds")

m = length(Kruskal_pvalue)
bonf_t = -log10(0.01 / m)
fst_cut = quantile(abs(FST_directional), 0.99, na.rm = TRUE)

FST_by_KW_test_no_zero = FST_by_KW_test[ !(FST_by_KW_test$Krusal_pvalue_directional == 0 & FST_by_KW_test$FST_directional == 0) , ]
figure = ggplot(FST_by_KW_test_no_zero, aes(x = Krusal_pvalue_directional, y = FST_directional, color = Color_the_figure)) +
  geom_point() + 
  theme_classic() +
  theme( legend.position = "bottom") +
  scale_color_viridis_c(direction = -1, option = "inferno") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_vline(xintercept =  c(-bonf_t, bonf_t),
             linetype   = "dotted",
             colour     = "brown") +
  geom_hline(yintercept =  c(-fst_cut, fst_cut),
             linetype   = "dotted", colour = "brown") +
  #xlab("-log10(KW p-value) SNV by AMY1 Copy") +
  xlab("") +
  #ylab("SNV FST") +
  ylab("") +
  xlim(-18, 18) +
  ylim(-0.42, 0.42) +
  labs(color = "FST by p-value")
  ggtitle("Alternative SNV Status: Chromosome 1")

pdf(file = "Cross_figure_Quechua_Maya.pdf", width = 10, height = 6)
print(figure)
dev.off()

