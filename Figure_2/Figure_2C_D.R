library(dplyr)
library(ggplot2)
library(furrr)
library(stringr)
plan(multisession, workers = 12)

Americas <- read.delim("Americas_all_genes.txt", header = TRUE)
Americas <- Americas %>% select(-Number)
tbl <- as_tibble(read.csv("tbl_protein_coding_genes.csv"))

# Filter populations and combine data
Seq1 <- Americas %>% dplyr::filter(Population == "QC_H")
Seq2 <- Americas %>% dplyr::filter(Population == "Maya")
test_combined_vst <- rbind(Seq1, Seq2)

# Subset the genes of interest (gene columns + Population)
shuffled_data = test_combined_vst

# Function to calculate Vst and combined metric
get_null_combined_metric <- function(shuffled_data) {
  Seq1 <- subset(shuffled_data, Population == "QC_H")
  Seq2 <- subset(shuffled_data, Population == "Maya")
  test_combined_vst <- rbind(Seq1, Seq2)
  
  # Function to calculate Vst
  getVst <- function(data) {
    if ((max(data) - min(data)) > 1) {
      dat1 <- data[1:nrow(Seq1)]
      dat2 <- data[(nrow(Seq1) + 1):(nrow(Seq1) + nrow(Seq2))]
      Vtotal <- var(c(dat1, dat2) )
      Vgroup <- ((var(dat1) * (length(dat1)-1)) +
                   (var(dat2) * (length(dat2)-1))) /
        ((length(dat1)-1) + (length(dat2)-1))
      Vst <- (Vtotal - Vgroup) / Vtotal
      if (Vtotal == 0 || Vst < 0) {
        Vst <- 0
      }
      return(Vst)
    } else {
      Vst <- 0
    }
  }
  
  # Extract gene column names
  gene_names <- colnames(test_combined_vst)[4:ncol(test_combined_vst)]
  
  # Calculate Vst for each gene
  Vst_values <- apply(test_combined_vst[4:ncol(test_combined_vst)], 2, getVst)
  
  # Combine Vst values and gene names into a data frame
  Vst <- data.frame(Gene = gene_names, VST = Vst_values)
  rownames(Vst) <- NULL
  
  # Calculate medians for each gene
  getMedians <- function(data) {
    dat1 <- data[1:nrow(Seq1)]
    dat2 <- data[(nrow(Seq1) + 1):(nrow(Seq1) + nrow(Seq2))]
    abs(median(dat1,) - median(dat2,))
  }
  medians <- apply(test_combined_vst[4:ncol(test_combined_vst)], 2, getMedians)
  
  Vst$Medians <- medians
  
  # Perform the inner join between Vst and tbl
  Vst <- Vst %>%
    inner_join(tbl, by = c("Gene" = "GENENAME"), relationship = "many-to-many") %>%
    group_by(Gene) %>%
    filter(VST == max(VST)) %>%
    ungroup() %>%
    distinct(Gene, .keep_all = TRUE)
  
  Outlier_removed <- Vst %>%
    dplyr::filter(
      Medians != 0 &
        VST > 0 &
        !str_detect(Gene, "USP17L") &
        !str_detect(Gene, "NBPF") &
        Gene != "MTRNR2L8" &
        !str_detect(Gene, "DUX4")
    )
  
  return(data.frame(
    Medians = Outlier_removed$Medians,
    VST = Outlier_removed$VST,
    Genes = Outlier_removed$Gene,
    stringsAsFactors = FALSE
  ))
}

set.seed(123)

# Generate null distribution
null_distribution <- future_map_dfr(1:10000, ~ {
  shuffled_temp <- shuffled_data
  shuffled_temp$Population <- sample(shuffled_temp$Population)
  
  # Calculate combined metrics
  result <- get_null_combined_metric(shuffled_temp)
  Medians <- result$Medians
  VST <- result$VST
  Genes <- result$Genes
  
  # Ensure result is a data frame before returning
  data.frame(Medians = Medians, VST = VST, Genes = Genes, stringsAsFactors = FALSE)
}, .options = furrr_options(seed = TRUE))

# Observed AMY1
observed_result <- get_null_combined_metric(test_combined_vst)
AMY1 = observed_result %>%
  dplyr::filter(Genes == "AMY1A")
print(paste("AMY1 Median:", AMY1$Medians))
print(paste("AMY1 Vst:", AMY1$VST))

null_distribution <- null_distribution %>%
  filter(!is.na(Medians), !is.na(VST))

# Calculate critical value (95th percentile) for the null distribution
critical_value_medians <- quantile(null_distribution$Medians, probs = 0.95)
critical_value_vst <- quantile(null_distribution$VST, probs = 0.95)

# Check if the observed value is significant
if (AMY1$Medians > critical_value_medians) {
  print("The top combined metric for medians is significant (p < 0.05).")
} else {
  print("The top combined metric for medians is not significant (p >= 0.05).")
}

if (AMY1$VST > critical_value_vst) {
  print("The top combined metric for VST is significant (p < 0.05).")
} else {
  print("The top combined metric for VST is not significant (p >= 0.05).")
}

# Visualization
QC_H_Maya_medians <- ggplot(data.frame(x = null_distribution$Medians), aes(x = x)) +
  geom_histogram(color = "#ffcc00ff", fill = "#ffcc0080") +
  xlab("All Medians 10000 Permutations") +
  ylab("Count") +
  geom_vline(xintercept = critical_value_medians, linewidth = 1.25, linetype = "dashed", color = "#cd1076ff") +
  geom_vline(xintercept = AMY1$Medians, linewidth = 1.25, linetype = "dashed", color = "#5d478bff") +
  theme_minimal()

QC_H_Maya_vst <- ggplot(data.frame(x = null_distribution$VST), aes(x = x)) +
  geom_histogram(color = "#ffcc00ff", fill = "#ffcc0080") +
  xlab("ALL VST values 10000 Permutations") +
  ylab("Count") +
  geom_vline(xintercept = critical_value_vst, linewidth = 1.25, linetype = "dashed", color = "#cd1076ff") +
  geom_vline(xintercept = AMY1$VST, linewidth = 1.25, linetype = "dashed", color = "#5d478bff") +
  theme_minimal()

# Save results
saveRDS(null_distribution, "All_median_vst_QC_H_Maya_5.rds")
pdf("QC_H_Maya_medians_10000.pdf", width = 6, height = 4)
print(QC_H_Maya_medians)
dev.off()

pdf("QC_H_Maya_VST_10000.pdf", width = 6, height = 4)
print(QC_H_Maya_vst)
dev.off()

