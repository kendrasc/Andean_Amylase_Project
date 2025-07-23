library(ggplot2)
library(dplyr)

population_names <- c("YRI","TSI", "STU", "PUR", "PJL", "PEL",
                      "MXL", "MSL", "LWK", "KHV", "JPT", "ITU",
                      "IBS", "GWD", "GIH", "GBR", "FIN", "ESN",
                      "CLM", "CHS", "CHB", "CEU", "CDX", "BEB",
                      "ASW", "ACB", "Quechua", "Maya_Chiapas")

setwd("~/Desktop/Amylase_Americas/SNP_pop_significance_comparison_figure")
GeneToCN <- read.delim("~/Desktop/Amylase_Americas/All_amylase_gene_copy_number.txt", header = TRUE)
results_df <- data.frame(Population = character(), Max_Pvalue = numeric(), SNP_Position = character(), stringsAsFactors = FALSE)

######################### Kruskal Wallis Analysis ##########################
# calculates KW for the SNPs --> sees what is the max p-value associated with AMY1 copy number for each pop
analyze_population = function(population_name) {
  
  if (population_name == "Quechua") {
    genotype_file = "Quechua_no_EUR_amylase_chr1_no_indels_0.05_no_missing_unique_Genotype.txt"
    GenetoCN_Pops = GeneToCN %>%
      filter(Pop %in% c("Highland_Quechua", "Lowland_Quechua")) %>%
      filter(!ID %in% c("Quechua_Lowland_3", "Quechua_Highland_23", "Quechua_Highland_24"))
  } else if (population_name == "Maya_Chiapas") {
    genotype_file = "Maya_Chiapas_amylase_chr1_no_indels_0.05_no_missing_unique_Genotype.txt"
    GenetoCN_Pops = GeneToCN %>%
      filter(Pop == population_name)
  } else {
    genotype_file = paste0("1k_genomes_", population_name, "_amylase_chr1_no_indels_0.05_no_missing_unique_Genotype.txt")
    GenetoCN_Pops = GeneToCN %>%
      filter(Pop == population_name)
  }
  
  Genotype = read.delim(genotype_file, header = TRUE) %>%
    select(-CHROM) %>%
    mutate(POS = paste0("X", POS))
  rownames(Genotype) = Genotype$POS
  
  Genotype_transposed = as.data.frame(t(Genotype[,-1])) %>%
    `colnames<-`(Genotype$POS) %>%
    mutate(ID = rownames(.)) %>%
    select(ID, everything())
  Genotype_transposed$ID = sub(".*_", "", Genotype_transposed$ID)
  
  GenetoCN_Pops_sorted = GenetoCN_Pops[order(GenetoCN_Pops$ID), ]
  Genotype_sorted = Genotype_transposed[order(Genotype_transposed$ID), ]
  dataset_for_snp_analyses = cbind(GenetoCN_Pops_sorted, Genotype_sorted[, -1])
  
  snp_data = dataset_for_snp_analyses[, 6:ncol(dataset_for_snp_analyses)]
  dataset_for_snp_analyses[6:ncol(dataset_for_snp_analyses)] <- lapply(
    snp_data,
    function(x) factor(
      ifelse(x == "0|0", "Homozygous Ancestral",
             ifelse(x == "1|1", "Homozygous Derived", "Heterozygous")),
      levels = c("Homozygous Ancestral", "Heterozygous", "Homozygous Derived")
    )
  )
  
  AMY1_data = dataset_for_snp_analyses$AMY1
  Kruskal_pvalue = vapply(snp_data, function(snp) {
    kruskal.test(AMY1_data ~ snp)$p.value
  }, numeric(1))
  
  max_pvalue = max(-log10(Kruskal_pvalue))
  max_snp_position = names(Kruskal_pvalue)[which.max(-log10(Kruskal_pvalue))]
  
  cat("Maximum p-value for", population_name, ":", max_pvalue, 
      "at", max_snp_position, "\n\n")
  results_df <<- rbind(results_df, data.frame(Population = population_name, Max_Pvalue = max_pvalue, SNP_Position = max_snp_position))
}
lapply(population_names, analyze_population)

###################### Results from LD analysis ###################
onek_genomes_LD = read.delim("~/Desktop/Amylase_Americas/Andean_LD_analysis/LD_file_103348464-103830994.txt")
# example here is 103348464-103830994 region

sorted_snps = results_df[order(results_df$Population),]
sorted_onek_genomes = onek_genomes_LD[order(onek_genomes_LD$Pop),]

sorted_onek_genomes$Max_Pvalue = sorted_snps$Max_Pvalue
sorted_onek_genomes$SNP_Position = sorted_snps$sorted_onek_genomes

# Example here = r^2 Across locus 
file = ggplot(sorted_onek_genomes, aes(x = Max_Pvalue, y = Mean_r2_Across, color = Region)) +
  geom_point(size = 4) +
  ylab("Average LD Across Chr1: 103348464-103830994") +
  xlab("Top Significantly Associated SNP with AMY1 Copy Number") +
  scale_color_manual(values = c("AFR" = "#cd1076bf", "EAS" = "#ffaaeebf", "AMR" = "#ffcc00bf",
                                "EUR" = "#ff6600bf", "SAS" = "#8b4726bf")) +
  theme_minimal()

pdf(file = "~/Desktop/Andean_Amylase_PDFs/Figure_S7.pdf", width = 8, height = 6)
print(file)
dev.off()

