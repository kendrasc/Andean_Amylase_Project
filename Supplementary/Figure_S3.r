library(dplyr)
library(ggplot2)
#library(ggrepel)
#library(AnnotationHub)
#library(ensembldb)
library(viridis)
library(stringr)

#################### All 1k genomes populations ########################
modern = read.csv("~/Desktop/modern_humans_redownload.csv")
n <- ncol(modern)
modern <- modern[, c(1, n-1, n, 2:(n-2))]


################## Protein coding genes #########################
tbl <- as_tibble(read.csv("~/Desktop/tbl_protein_coding_genes.csv"))
############# Pulls all genes from Human Ensembl #######################
#hub = AnnotationHub()
#query(hub, c("EnsDb", "Homo sapiens", "97"))
#edb = hub[["AH73881"]]
#keytypes(edb)
#columns(edb)
#keys = keys(edb, "GENENAME")
#columns =  c("GENEID", "ENTREZID", "GENEBIOTYPE")
#tbl =
#  ensembldb::select(edb, keys, columns, keytype = "GENENAME") %>%
#  as_tibble()

################ This is filtering out the amy genes that have bigger values = fix here
######### Pull Only Protein Coding Genes (23069) #######################
#supportedFilters()
#filter = ~ gene_name %in% keys & gene_biotype == "protein_coding"
#tbl =
#  ensembldb::select(edb, filter, columns) %>%
#  as_tibble()

#write.csv(tbl, "~/Desktop/tbl_protein_coding_genes.csv", row.names = FALSE)


######################################### Calculating VST ################################################################################################
Seq1 = modern %>% dplyr::filter(pop == "PEL")
Seq2 = modern %>% dplyr::filter(pop == "MXL")
test_combined_vst = rbind(Seq1, Seq2)

getVst <- function(data) {
  if ((max(data) - min(data)) > 1){
    dat1 <- data[1:(length(Seq1$pop))]
    dat2 <- data[(length(Seq1$pop)+1):(length(Seq1$pop) + (length(Seq2$pop)))]
    Vtotal <- var(c(dat1,dat2))
    Vgroup <- ((var(dat1) * (length(dat1)-1)) + (var(dat2) * (length(dat2)-1))) /
      ((length(dat1)-1) + (length(dat2)-1))
    Vst <- c((Vtotal-Vgroup) / Vtotal)
    if (Vst < 0) {
      Vst = 0
    }
    return(Vst)
  }
  else{
    Vst = 0
  }
}

full_dataset = as.data.frame(apply(test_combined_vst[4:ncol(test_combined_vst)],2, getVst))
full_dataset$Gene = row.names(full_dataset)
rownames(full_dataset) <- NULL
colnames(full_dataset) = c("VST", "Gene")

######################################### Calculate Medians #############################
Seq1_median = apply(Seq1[4:ncol(test_combined_vst)], 2, median)
Seq2_median = apply(Seq2[4:ncol(test_combined_vst)], 2, median)
median_across = as.data.frame(abs(Seq1_median - Seq2_median))
median_across$Gene = row.names(median_across)
rownames(median_across) <- NULL
colnames(median_across) = c("Medians", "Gene")

full_dataset$Medians = median_across$Medians
######################################### Filter only the Top VST Hits ##############################
full_dataset$Gene <- sub("\\..*", "", full_dataset$Gene)

# Perform the inner join between full_dataset and tbl based on Gene and GENENAME
filtered_dataset <- full_dataset %>%
  inner_join(tbl, by = c("Gene" = "GENENAME"), relationship = "many-to-many") %>%
  group_by(Gene) %>%
  # For each Gene, keep only the row with the largest VST value
  dplyr::filter(VST == max(VST)) %>%
  ungroup() %>%
  distinct(Gene, .keep_all = TRUE)

######################################### VST Figure #####################################################################################################
locations = read.delim("~/Desktop/Stuff_on_Desktop/Gene_CNVs/locations_please_work.txt", header = T)

# file to get unique locations --> Choosing the first location
locations_unique = locations %>%
  dplyr::distinct(across(Gene),
                  .keep_all = TRUE)

locations_adjusted = locations_unique[locations_unique$Gene %in% filtered_dataset$Gene,]
locations_adjusted_sorted <- locations_adjusted[order(locations_adjusted$Gene),]
filtered_dataset_sorted <- filtered_dataset[order(filtered_dataset$Gene),]

############################### Filter out usually high copy numbered genes/clear errors #############
Outlier_removed <- filtered_dataset_sorted %>%
  dplyr::filter(
    Medians != 0 &
      VST != 0 &
      !str_detect(Gene, "USP17L") &
      !str_detect(Gene, "NBPF") &
      Gene != "MTRNR2L8" &
      !str_detect(Gene, "DUX4")
  )

Outlier_removed$combined_metric <- Outlier_removed$Medians * Outlier_removed$VST

PEL = ggplot(Outlier_removed, aes(x = Medians, y = VST, color = combined_metric)) +
  geom_point(size = 3) +
  scale_color_gradientn(
    colors = rev(viridis(256, option = "inferno")), # Reverses the color order
    limits = c(0, 1),
    oob = scales::squish) +
  ylim(0, 0.35) +
  xlim(0, 11) +
  xlab("Median Gene Copy Number Difference") + 
  ylab("VST") + 
  ggtitle("PEL vs. MXL") + 
  labs(color = " VST x Median") +
  theme_minimal() +
  theme(legend.position = "bottom")

write.csv(Outlier_removed, "~/Desktop/Amylase_Americas/CSVs/PEL_vs_MXL_VST_positive_VST_values_fixed_medians.csv")

pdf(file = "~/Desktop/Amylase_Americas/PDFs/PEL_MXL_median_VST.pdf", width = 8, height = 4.5)
print(PEL)
dev.off()
