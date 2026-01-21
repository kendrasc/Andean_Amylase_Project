library(tidyverse)
library(viridis)
library(ggpubr)
library(ggbreak)

setwd("~/Desktop/Amylase_Americas/")
GeneToCN = read.delim("All_amylase_gene_copy_number.txt", header = T)
Unrelated = read.delim("Unrelated_labeled_1k_genomes.txt", header = F) # Just an extra check
split_vals = do.call(rbind, strsplit(Unrelated$V1, " "))
Unrelated = data.frame(ID = split_vals[, 1], Pop = split_vals[, 2])

pop_list = c("CDX", "CHB", "KHV", "JPT", "CHS",
              "CEU", "GBR", "TSI", "FIN", "IBS", 
              "BEB", "GIH", "PJL", "STU","ITU", 
              "YRI", "LWK", "ESN", "MSL", "GWD", "ACB", "ASW",
              "PEL", "MXL", "PUR", "CLM" )

global_counts = c()
plot_list = list()
df_list = list()
dunn_list = list()

for (pop_name in pop_list) {
  message("Processing: ", pop_name)
  
Genotype = read.delim(
  paste0("1k_genomes_", pop_name, "_amylase_chr1_no_indels_0.05_no_missing_unique_Genotype.txt"),
  header = TRUE)

Genotype %>%
  dplyr::filter(POS > 103348464 & POS < 103830994)

Genotype = subset(Genotype, select = -CHROM)
Genotype = Genotype %>%
  mutate(POS = paste0("X", POS))
rownames(Genotype) = Genotype$POS
Genotype_transposed = as.data.frame(t(Genotype[,-1])) 
colnames(Genotype_transposed) = Genotype$POS
Genotype_transposed$ID = rownames(Genotype_transposed) 
rownames(Genotype_transposed) = NULL
Genotype_transposed = Genotype_transposed %>%
  select(ID, everything()) %>%
  dplyr::filter(ID %in% Unrelated$ID)

GenetoCN_Pops = GeneToCN %>%
  dplyr::filter(Pop == pop_name) %>%
  dplyr::filter(ID %in% Unrelated$ID)

################################# SNP by Gene Copy Number ###############################
GenetoCN_Pops_sorted = GenetoCN_Pops[order(GenetoCN_Pops$ID),]
Genotype_sorted = Genotype_transposed[order(Genotype_transposed$ID),]
dataset_for_snp_analyses = cbind(GenetoCN_Pops_sorted, Genotype_sorted[, -1])

# Relabel SNP genotypes and save
dataset_for_snp_analyses[8:ncol(dataset_for_snp_analyses)] = lapply(
  dataset_for_snp_analyses[8:ncol(dataset_for_snp_analyses)],
  function(x) factor(ifelse(x == "0|0", "Homozygous Ancestral",
                            ifelse(x == "1|1", "Homozygous Derived", "Heterozygous")),
                     levels = c("Homozygous Ancestral", "Heterozygous", "Homozygous Derived"))
)

dataset_for_snp_analyses$AMY1 = round(dataset_for_snp_analyses$AMY1)

# Apply the Kruskal-Wallis test across the columns and extract p-values 
AMY1_data = dataset_for_snp_analyses$AMY1
snp_data = dataset_for_snp_analyses[, 8:length(dataset_for_snp_analyses)]

Kruskal_pvalue = vapply(snp_data, function(snp) {
  kruskal.test(AMY1_data ~ snp)$p.value
}, numeric(1))

KW_adj = p.adjust(Kruskal_pvalue, method = "bonferroni") 

sig_snps_names = names(KW_adj)[KW_adj < 0.05]
sig_snps = KW_adj[KW_adj < 0.05]

########### Giving a direction ########

Kruskal_Wallis_status = numeric()
Kruskal_Wallis_status = sapply(dataset_for_snp_analyses[, 9:length(dataset_for_snp_analyses)], function(col) {
  
  # Filter values in gene copy number column based on current column
  Hets = dataset_for_snp_analyses[[3]][col == "Heterozygous"]
  Ancestral = dataset_for_snp_analyses[[3]][col == "Homozygous Ancestral"]
  Derived = dataset_for_snp_analyses[[3]][col == "Homozygous Derived"]
  
  # Calculate medians for the filtered groups, handling empty groups
  Het_median = if (length(Hets) > 0) median(Hets) else NA
  Ancestral_median = if (length(Ancestral) > 0) median(Ancestral) else NA
  Derived_median = if (length(Derived) > 0) median(Derived) else NA
  
  # Set conditions based on medians, ensuring no NA values in comparisons
  if (!is.na(Derived_median) && !is.na(Het_median) && !is.na(Ancestral_median)) {
    condition_pos_met = Derived_median > Het_median && Het_median > Ancestral_median
    condition_neg_met = Derived_median < Het_median && Het_median < Ancestral_median
  } else if (!is.na(Derived_median) && !is.na(Het_median) && is.na(Ancestral_median)) {
    condition_pos_met = Derived_median > Het_median
    condition_neg_met = Derived_median < Het_median
  } else if (is.na(Derived_median) && !is.na(Het_median) && !is.na(Ancestral_median)) {
    condition_pos_met = Het_median > Ancestral_median
    condition_neg_met = Het_median < Ancestral_median
  } else if (!is.na(Derived_median) && is.na(Het_median) && !is.na(Ancestral_median)) {
    condition_pos_met = Derived_median > Ancestral_median
    condition_neg_met = Derived_median < Ancestral_median
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

Kruskal_Wallis_status_adj = as.numeric(unlist(Kruskal_Wallis_status))
Kruskal_Wallis_status_SNV_IDs = names(Kruskal_Wallis_status)
Kruskal_Wallis_status_adj_named = data.frame(Kruskal_Wallis_status_SNV_IDs, Kruskal_Wallis_status_adj)
Kruskal_Wallis_status_sig_snps = Kruskal_Wallis_status_adj_named %>%
  dplyr::filter(Kruskal_Wallis_status_SNV_IDs %in% sig_snps_names)


########### Calculating Frequencies
sig_snps_genotypes = dataset_for_snp_analyses %>% select(all_of(sig_snps_names))

Frequencies = sapply(sig_snps_names, function(col_name) {
  
  col = dataset_for_snp_analyses[[col_name]]
  count = table(factor(col, levels = c("Heterozygous", "Homozygous Ancestral", "Homozygous Derived")))
  
  # Extract counts
  Hets = unname(count["Heterozygous"])
  Ancestral = unname(count["Homozygous Ancestral"])
  Derived = unname(count["Homozygous Derived"])
  
  # Replace NAs with 0
  Hets = ifelse(is.na(Hets), 0, Hets)
  Ancestral = ifelse(is.na(Ancestral), 0, Ancestral)
  Derived = ifelse(is.na(Derived), 0, Derived)
  
  # Calculate derived allele frequency
  total = Hets + Ancestral + Derived
  if (total == 0) {
    return(0)
  }
  Frequency_derived = (Hets + 2 * Derived) / (2 * total)
  return(Frequency_derived)
})
  
############# Main figure
Krusal_pvalue_bonferroni = Kruskal_Wallis_status_sig_snps$Kruskal_Wallis_status_adj*(-log10(as.numeric(unlist(sig_snps))))
Frequency = as.numeric(unlist(Frequencies)) # frequency of derived state
SNV_IDs = names(Frequencies)
# gives me IDs, the KW values + which direction association + freq of derived state
data_frame_version = data.frame(SNV_IDs, Krusal_pvalue_bonferroni, Frequency)

idx = data_frame_version$Frequency < 0.5

data_frame_version$Frequency[idx] = -1*(1 - data_frame_version$Frequency[idx])
data_frame_version$Krusal_pvalue_bonferroni[idx] = -data_frame_version$Krusal_pvalue_bonferroni[idx] 

data_frame_version = data_frame_version %>%
  mutate(Pop = pop_name)
df_list[[pop_name]] = data_frame_version 

label_positive_derived = paste("Count:", sum(data_frame_version$Krusal_pvalue_bonferroni > 0 & data_frame_version$Frequency > 0))
label_negative_derived = paste("Count:", sum(data_frame_version$Krusal_pvalue_bonferroni < 0 & data_frame_version$Frequency > 0))
label_positive_ancestral = paste("Count:", sum(data_frame_version$Krusal_pvalue_bonferroni > 0 & data_frame_version$Frequency < 0))
label_negative_ancestral = paste("Count:", sum(data_frame_version$Krusal_pvalue_bonferroni < 0 & data_frame_version$Frequency < 0))


plot_data = data_frame_version %>% 
  mutate(
    side = ifelse(Frequency < 0, "reference", "alternate"),
    side = factor(side, levels = c("reference", "alternate"))) %>% 
  count(side, Frequency, Krusal_pvalue_bonferroni, name = "n")

global_counts = c(global_counts, plot_data$n)

label_df = tibble::tribble(
  ~side,        ~Frequency, ~Krusal_pvalue_bonferroni,  ~label,                   ~hjust,
  "alternate",   1,          7.5,  label_positive_derived,   1,
  "reference",  -1,          7.5,  label_positive_ancestral, 0,
  "alternate",   1,         -7.5,  label_negative_derived,   1,
  "reference",  -1,         -7.5,  label_negative_ancestral, 0) %>% 
  mutate(side = factor(side, levels = c("reference", "alternate")))

boundary_df = tibble::tribble(
  ~side,       ~Frequency, ~Krusal_pvalue_bonferroni,
  "reference", -1.0,        0,
  "reference", -0.5,        0,
  "alternate",  0.5,        0,
  "alternate",  1.0,        0
) |>
  mutate(side = factor(side, levels = c("reference", "alternate")))

figure = ggplot(plot_data,
               aes(x = Frequency,
                   y = Krusal_pvalue_bonferroni,
                   colour = Krusal_pvalue_bonferroni,
                   size   = n,
                   alpha  = 0.75)) +
  geom_blank(data = boundary_df, inherit.aes = FALSE,
             aes(x = Frequency, y = Krusal_pvalue_bonferroni)) +
  geom_point() +
  geom_text(data = label_df,
            aes(x = Frequency, y = Krusal_pvalue_bonferroni,
                label = label, hjust = hjust),
            inherit.aes = FALSE, vjust = 1, size = 4) +
  scale_color_viridis(limits = c(-8, 8), option = "inferno") +
  scale_size_continuous(limits = global_counts, range = c(1, 6)) +
  scale_y_continuous(limits = c(-8, 8)) +
  facet_grid(. ~ side, scales = "free_x", space = "free_x") +
  theme_minimal() +
  theme(axis.title   = element_blank(),
        panel.spacing.x = unit(1, "lines"),
        panel.clip = "off") +
  ggtitle(pop_name)

plot_list[[pop_name]] = figure
}

#################### Quechua ################

pop_name = "Quechua"

Genotype = read.delim("Quechua_PhasedChr1_no_missing_biallelic_no_EUR_maf_0.05_genotypes_recreate_Genotype.txt", header = T)

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

GenetoCN_Pops = GeneToCN %>%
  dplyr::filter(Pop %in% c("Lowland_Andeans", "Highland_Andeans")) %>%
# remove three individuals with non-American Ancestry here
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

# Relabel SNP genotypes and save
dataset_for_snp_analyses[9:ncol(dataset_for_snp_analyses)] = lapply(
  dataset_for_snp_analyses[9:ncol(dataset_for_snp_analyses)],
  function(x) factor(ifelse(x == "0|0", "Homozygous Ancestral",
                            ifelse(x == "1|1", "Homozygous Derived", "Heterozygous")),
                     levels = c("Homozygous Ancestral", "Heterozygous", "Homozygous Derived")))


dataset_for_snp_analyses$AMY1 = round(dataset_for_snp_analyses$AMY1)

# Apply the Kruskal-Wallis test across the columns and extract p-values 
AMY1_data = dataset_for_snp_analyses$AMY1
snp_data = dataset_for_snp_analyses[, 9:length(dataset_for_snp_analyses)]

Kruskal_pvalue = vapply(snp_data, function(snp) {
  kruskal.test(AMY1_data ~ snp)$p.value
}, numeric(1))

KW_check = as.numeric(unlist(Kruskal_pvalue))
SNV_IDs_check = names(Kruskal_pvalue)


KW_adj = p.adjust(Kruskal_pvalue, method = "bonferroni") 

sig_snps_names = names(KW_adj)[KW_adj < 0.05]
sig_snps = KW_adj[KW_adj < 0.05]

########### Giving a direction ########

Kruskal_Wallis_status = numeric()
Kruskal_Wallis_status = sapply(dataset_for_snp_analyses[, 9:length(dataset_for_snp_analyses)], function(col) {
  
  # Filter values in gene copy number column based on current column
  Hets = dataset_for_snp_analyses[[3]][col == "Heterozygous"]
  Ancestral = dataset_for_snp_analyses[[3]][col == "Homozygous Ancestral"]
  Derived = dataset_for_snp_analyses[[3]][col == "Homozygous Derived"]
  
  # Calculate medians for the filtered groups, handling empty groups
  Het_median = if (length(Hets) > 0) median(Hets) else NA
  Ancestral_median = if (length(Ancestral) > 0) median(Ancestral) else NA
  Derived_median = if (length(Derived) > 0) median(Derived) else NA
  
  # DeMayae conditions based on medians, ensuring no NA values in comparisons
  if (!is.na(Derived_median) && !is.na(Het_median) && !is.na(Ancestral_median)) {
    condition_pos_met = Derived_median > Het_median && Het_median > Ancestral_median
    condition_neg_met = Derived_median < Het_median && Het_median < Ancestral_median
  } else if (!is.na(Derived_median) && !is.na(Het_median) && is.na(Ancestral_median)) {
    condition_pos_met = Derived_median > Het_median
    condition_neg_met = Derived_median < Het_median
  } else if (is.na(Derived_median) && !is.na(Het_median) && !is.na(Ancestral_median)) {
    condition_pos_met = Het_median > Ancestral_median
    condition_neg_met = Het_median < Ancestral_median
  } else if (!is.na(Derived_median) && is.na(Het_median) && !is.na(Ancestral_median)) {
    condition_pos_met = Derived_median > Ancestral_median
    condition_neg_met = Derived_median < Ancestral_median
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

Kruskal_Wallis_status_adj = as.numeric(unlist(Kruskal_Wallis_status))
Kruskal_Wallis_status_SNV_IDs = names(Kruskal_Wallis_status)
Kruskal_Wallis_status_adj_named = data.frame(Kruskal_Wallis_status_SNV_IDs, Kruskal_Wallis_status_adj)
Kruskal_Wallis_status_sig_snps = Kruskal_Wallis_status_adj_named %>%
  dplyr::filter(Kruskal_Wallis_status_SNV_IDs %in% sig_snps_names)


########### Calculating Frequencies
sig_snps_genotypes = dataset_for_snp_analyses %>% select(all_of(sig_snps_names))

Frequencies = sapply(sig_snps_names, function(col_name) {
  
  col = dataset_for_snp_analyses[[col_name]]
  count = table(factor(col, levels = c("Heterozygous", "Homozygous Ancestral", "Homozygous Derived")))
  
  # Extract counts
  Hets = unname(count["Heterozygous"])
  Ancestral = unname(count["Homozygous Ancestral"])
  Derived = unname(count["Homozygous Derived"])
  
  # Replace NAs with 0
  Hets = ifelse(is.na(Hets), 0, Hets)
  Ancestral = ifelse(is.na(Ancestral), 0, Ancestral)
  Derived = ifelse(is.na(Derived), 0, Derived)
  
  # Calculate derived allele frequency
  total = Hets + Ancestral + Derived
  if (total == 0) {
    return(0)
  }
  Frequency_derived = (Hets + 2 * Derived) / (2 * total)
  return(Frequency_derived)
})

############################################

dunn_values = lapply(sig_snps_names, function(snp) {
  res = dunnTest(AMY1_data ~ dataset_for_snp_analyses[[snp]],
                  method = "bonferroni")$res
  cbind(SNP = snp,
        KW_q = KW_adj[[snp]],
        res)
})

dunn_results = do.call(rbind, dunn_values)
if (!is.null(dunn_results) && nrow(dunn_results) > 0) {
  dunn_results$Pop = pop_name
  dunn_list[[pop_name]] <- dunn_results
}


###################### Figure
Krusal_pvalue_bonferroni = Kruskal_Wallis_status_sig_snps$Kruskal_Wallis_status_adj*(-log10(as.numeric(unlist(sig_snps))))
Frequency = as.numeric(unlist(Frequencies)) # frequency of derived state
SNV_IDs = names(Frequencies)
# gives me IDs, the KW values + which direction association + freq of derived state
data_frame_version = data.frame(SNV_IDs, Krusal_pvalue_bonferroni, Frequency)

idx = data_frame_version$Frequency < 0.5
data_frame_version$Frequency[idx] = -1*(1 - data_frame_version$Frequency[idx])
data_frame_version$Krusal_pvalue_bonferroni[idx] = -data_frame_version$Krusal_pvalue_bonferroni[idx] 

data_frame_version = data_frame_version %>%
  mutate(Pop = pop_name)
df_list[[pop_name]] = data_frame_version 

label_positive_derived = paste("Count:", sum(data_frame_version$Krusal_pvalue_bonferroni > 0 & data_frame_version$Frequency > 0))
label_negative_derived = paste("Count:", sum(data_frame_version$Krusal_pvalue_bonferroni < 0 & data_frame_version$Frequency > 0))
label_positive_ancestral = paste("Count:", sum(data_frame_version$Krusal_pvalue_bonferroni > 0 & data_frame_version$Frequency < 0))
label_negative_ancestral = paste("Count:", sum(data_frame_version$Krusal_pvalue_bonferroni < 0 & data_frame_version$Frequency < 0))

plot_data = data_frame_version %>% 
  mutate(
    side = ifelse(Frequency < 0, "reference", "alternate"),
    side = factor(side, levels = c("reference", "alternate"))
  ) %>% 
  count(side, Frequency, Krusal_pvalue_bonferroni, name = "n")

global_counts = c(global_counts, plot_data$n)

label_df = tibble::tribble(
  ~side,        ~Frequency, ~Krusal_pvalue_bonferroni,  ~label,                   ~hjust,
  "alternate",   1,          7.5,  label_positive_derived,   1,
  "reference",  -1,          7.5,  label_positive_ancestral, 0,
  "alternate",   1,         -7.5,  label_negative_derived,   1,
  "reference",  -1,         -7.5,  label_negative_ancestral, 0
) %>% 
  mutate(side = factor(side, levels = c("reference", "alternate")))

boundary_df = tibble::tribble(
  ~side,       ~Frequency, ~Krusal_pvalue_bonferroni,
  "reference", -1.0,        0,
  "reference", -0.5,        0,
  "alternate",  0.5,        0,
  "alternate",  1.0,        0
) |>
  mutate(side = factor(side, levels = c("reference", "alternate")))

Quechua = ggplot(plot_data,
               aes(x = Frequency,
                   y = Krusal_pvalue_bonferroni,
                   colour = Krusal_pvalue_bonferroni,
                   size   = n,
                   alpha  = 0.75)) +
  geom_blank(data = boundary_df, inherit.aes = FALSE,
             aes(x = Frequency, y = Krusal_pvalue_bonferroni)) +
  geom_point() +
  geom_text(data = label_df,
            aes(x = Frequency, y = Krusal_pvalue_bonferroni,
                label = label, hjust = hjust),
            inherit.aes = FALSE, vjust = 1, size = 4) +
  scale_color_viridis(limits = c(-8, 8), option = "inferno") +
  scale_size_continuous(limits = global_counts, range = c(1, 6)) +
  scale_y_continuous(limits = c(-8, 8)) +
  facet_grid(. ~ side, scales = "free_x", space = "free_x") +
  theme_minimal() +
  theme(axis.title   = element_blank(),
        panel.spacing.x = unit(1, "lines"),
        panel.clip = "off") +
  ggtitle(pop_name)

####################### Maya ###################
pop_name = "Maya"
Genotype = read.delim("Maya_PhasedChr1_no_missing_biallelic_no_EUR_maf_0.05_genotypes_recreate_Genotype.txt", header = T)

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

GenetoCN_Pops = GeneToCN %>%
  dplyr::filter(Pop %in% c("Maya_Chiapas"))

################################# SNP by Gene Copy Number ###############################
GenetoCN_Pops_sorted = GenetoCN_Pops[order(GenetoCN_Pops$ID),]

Genotype_transposed = Genotype_transposed %>% 
  mutate(ID_short = sub(".*_", "", ID)) %>% 
  select(ID_short, everything())

Genotype_sorted = Genotype_transposed[order(Genotype_transposed$ID_short),]
dataset_for_snp_analyses = cbind(GenetoCN_Pops_sorted, Genotype_sorted[, -1])

# Relabel SNP genotypes and save
dataset_for_snp_analyses[9:ncol(dataset_for_snp_analyses)] = lapply(
  dataset_for_snp_analyses[9:ncol(dataset_for_snp_analyses)],
  function(x) factor(ifelse(x == "0|0", "Homozygous Ancestral",
                            ifelse(x == "1|1", "Homozygous Derived", "Heterozygous")),
                     levels = c("Homozygous Ancestral", "Heterozygous", "Homozygous Derived"))
)

dataset_for_snp_analyses$AMY1 = round(dataset_for_snp_analyses$AMY1)

# Apply the Kruskal-Wallis test across the columns and extract p-values 
AMY1_data = dataset_for_snp_analyses$AMY1
snp_data = dataset_for_snp_analyses[, 9:length(dataset_for_snp_analyses)]

Kruskal_pvalue = vapply(snp_data, function(snp) {
  kruskal.test(AMY1_data ~ snp)$p.value
}, numeric(1))

KW_adj = p.adjust(Kruskal_pvalue, method = "bonferroni") 

sig_snps_names = names(KW_adj)[KW_adj < 0.05]
sig_snps = KW_adj[KW_adj < 0.05]

########### Giving a direction ########

Kruskal_Wallis_status = numeric()
Kruskal_Wallis_status = sapply(dataset_for_snp_analyses[, 9:length(dataset_for_snp_analyses)], function(col) {
  
  # Filter values in gene copy number column based on current column
  Hets = dataset_for_snp_analyses[[3]][col == "Heterozygous"]
  Ancestral = dataset_for_snp_analyses[[3]][col == "Homozygous Ancestral"]
  Derived = dataset_for_snp_analyses[[3]][col == "Homozygous Derived"]
  
  # Calculate medians for the filtered groups, handling empty groups
  Het_median = if (length(Hets) > 0) median(Hets) else NA
  Ancestral_median = if (length(Ancestral) > 0) median(Ancestral) else NA
  Derived_median = if (length(Derived) > 0) median(Derived) else NA
  
  # DeMayae conditions based on medians, ensuring no NA values in comparisons
  if (!is.na(Derived_median) && !is.na(Het_median) && !is.na(Ancestral_median)) {
    condition_pos_met = Derived_median > Het_median && Het_median > Ancestral_median
    condition_neg_met = Derived_median < Het_median && Het_median < Ancestral_median
  } else if (!is.na(Derived_median) && !is.na(Het_median) && is.na(Ancestral_median)) {
    condition_pos_met = Derived_median > Het_median
    condition_neg_met = Derived_median < Het_median
  } else if (is.na(Derived_median) && !is.na(Het_median) && !is.na(Ancestral_median)) {
    condition_pos_met = Het_median > Ancestral_median
    condition_neg_met = Het_median < Ancestral_median
  } else if (!is.na(Derived_median) && is.na(Het_median) && !is.na(Ancestral_median)) {
    condition_pos_met = Derived_median > Ancestral_median
    condition_neg_met = Derived_median < Ancestral_median
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

Kruskal_Wallis_status_adj = as.numeric(unlist(Kruskal_Wallis_status))
Kruskal_Wallis_status_SNV_IDs = names(Kruskal_Wallis_status)
Kruskal_Wallis_status_adj_named = data.frame(Kruskal_Wallis_status_SNV_IDs, Kruskal_Wallis_status_adj)
Kruskal_Wallis_status_sig_snps = Kruskal_Wallis_status_adj_named %>%
  dplyr::filter(Kruskal_Wallis_status_SNV_IDs %in% sig_snps_names)


########### Calculating Frequencies

sig_snps_genotypes = dataset_for_snp_analyses %>% select(all_of(sig_snps_names))

Frequencies = sapply(sig_snps_names, function(col_name) {
  
  col = dataset_for_snp_analyses[[col_name]]
  count = table(factor(col, levels = c("Heterozygous", "Homozygous Ancestral", "Homozygous Derived")))
  
  # Extract counts
  Hets = unname(count["Heterozygous"])
  Ancestral = unname(count["Homozygous Ancestral"])
  Derived = unname(count["Homozygous Derived"])
  
  # Replace NAs with 0
  Hets = ifelse(is.na(Hets), 0, Hets)
  Ancestral = ifelse(is.na(Ancestral), 0, Ancestral)
  Derived = ifelse(is.na(Derived), 0, Derived)
  
  # Calculate derived allele frequency
  total = Hets + Ancestral + Derived
  if (total == 0) {
    return(0)
  }
  Frequency_derived = (Hets + 2 * Derived) / (2 * total)
  return(Frequency_derived)
})

############################################

dunn_values = lapply(sig_snps_names, function(snp) {
  res = dunnTest(AMY1_data ~ dataset_for_snp_analyses[[snp]],
                 method = "bonferroni")$res
  cbind(SNP = snp,
        KW_q = KW_adj[[snp]],
        res)
})

dunn_results = do.call(rbind, dunn_values)
if (!is.null(dunn_results) && nrow(dunn_results) > 0) {
  dunn_results$Pop = pop_name
  dunn_list[[pop_name]] <- dunn_results
}


###################### Figure
Krusal_pvalue_bonferroni = Kruskal_Wallis_status_sig_snps$Kruskal_Wallis_status_adj*(-log10(as.numeric(unlist(sig_snps))))
Frequency = as.numeric(unlist(Frequencies)) # frequency of derived state
SNV_IDs = names(Frequencies)
# gives me IDs, the KW values + which direction association + freq of derived state
data_frame_version = data.frame(SNV_IDs, Krusal_pvalue_bonferroni, Frequency)

idx = data_frame_version$Frequency < 0.5
data_frame_version$Frequency[idx] = -1*(1 - data_frame_version$Frequency[idx])
data_frame_version$Krusal_pvalue_bonferroni[idx] = -data_frame_version$Krusal_pvalue_bonferroni[idx] 

data_frame_version = data_frame_version %>%
  mutate(Pop = pop_name)
# add a Pop column
df_list[[pop_name]] = data_frame_version 

label_positive_derived = paste("Count:", sum(data_frame_version$Krusal_pvalue_bonferroni > 0 & data_frame_version$Frequency > 0))
label_negative_derived = paste("Count:", sum(data_frame_version$Krusal_pvalue_bonferroni < 0 & data_frame_version$Frequency > 0))
label_positive_ancestral = paste("Count:", sum(data_frame_version$Krusal_pvalue_bonferroni > 0 & data_frame_version$Frequency < 0))
label_negative_ancestral = paste("Count:", sum(data_frame_version$Krusal_pvalue_bonferroni < 0 & data_frame_version$Frequency < 0))

plot_data = data_frame_version %>% 
  mutate(
    side = ifelse(Frequency < 0, "reference", "alternate"),
    side = factor(side, levels = c("reference", "alternate"))
  ) %>% 
  count(side, Frequency, Krusal_pvalue_bonferroni, name = "n")

global_counts = c(global_counts, plot_data$n)

label_df = tibble::tribble(
  ~side,        ~Frequency, ~Krusal_pvalue_bonferroni,  ~label,                   ~hjust,
  "alternate",   1,          7.5,  label_positive_derived,   1,
  "reference",  -1,          7.5,  label_positive_ancestral, 0,
  "alternate",   1,         -7.5,  label_negative_derived,   1,
  "reference",  -1,         -7.5,  label_negative_ancestral, 0
) %>% 
  mutate(side = factor(side, levels = c("reference", "alternate")))


boundary_df = tibble::tribble(
  ~side,       ~Frequency, ~Krusal_pvalue_bonferroni,
  "reference", -1.0,        0,
  "reference", -0.5,        0,
  "alternate",  0.5,        0,
  "alternate",  1.0,        0
) |>
  mutate(side = factor(side, levels = c("reference", "alternate")))

Maya = ggplot(plot_data,
                 aes(x = Frequency,
                     y = Krusal_pvalue_bonferroni,
                     colour = Krusal_pvalue_bonferroni,
                     size   = n,
                     alpha  = 0.75)) +
  geom_blank(data = boundary_df, inherit.aes = FALSE,
             aes(x = Frequency, y = Krusal_pvalue_bonferroni)) +
  geom_point() +
  geom_text(data = label_df,
            aes(x = Frequency, y = Krusal_pvalue_bonferroni,
                label = label, hjust = hjust),
            inherit.aes = FALSE, vjust = 1, size = 4) +
  scale_color_viridis(limits = c(-8, 8), option = "inferno") +
  scale_size_continuous(limits = global_counts, range = c(1, 6)) +
  scale_y_continuous(limits = c(-8, 8)) +
  facet_grid(. ~ side, scales = "free_x", space = "free_x") +
  theme_minimal() +
  theme(axis.title   = element_blank(),
        panel.spacing.x = unit(1, "lines"),
        panel.clip = "off") +
  ggtitle(pop_name)


################# Make finalized figure ################

global_size_limits = range(global_counts, na.rm = TRUE)

plot_list[["Quechua"]] = Quechua
plot_list[["Maya"]] = Maya

plot_list = lapply(
  plot_list,
  `+`,
  scale_size_continuous(limits = global_size_limits, range = c(1, 6))
)

combined_plot = ggarrange(
  plotlist = plot_list,
  ncol     = 5,
  nrow     = ceiling(length(plot_list) / 5),
  common.legend = TRUE)

combined_plot_labeled = annotate_figure(
  combined_plot,
  left = text_grob("-log10 KW Bonferroni Corrected p-value by AMY1 Copy Association", rot = 90, size = 20),
  bottom = text_grob("Frequency of SNV in Population", size = 20))

ggsave("Combined_KW_histograms_by_pop.pdf", combined_plot_labeled, width = 18, height = 18)

big_df = dplyr::bind_rows(df_list)

write.csv(big_df, "Significant_SNVs_with_AMY1.csv")

dunn_output = dplyr::bind_rows(dunn_list)
write.csv(dunn_output, "Dunn_test_results.csv")
