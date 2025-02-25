library(dplyr)
library(ggplot2)
library(AnnotationHub)
library(ensembldb)
library(viridis)
library(stringr)

################# American  genomes populations ########################
Americas = read.delim("~/read/mrcanavar/data/Americas_all_genes.txt", header = TRUE, sep = "\t")

############# Pulls all genes from Human Ensembl #######################
hub = AnnotationHub()
query(hub, c("EnsDb", "Homo sapiens", "97"))
edb = hub[["AH73881"]]
keytypes(edb)
columns(edb)
keys = keys(edb, "GENENAME")
columns =  c("GENEID", "ENTREZID", "GENEBIOTYPE")
tbl =
  ensembldb::select(edb, keys, columns, keytype = "GENENAME") %>%
  as_tibble()

######### Pull Only Protein Coding Genes ###############
supportedFilters()
filter = ~ gene_name %in% keys & gene_biotype == "protein_coding"
tbl =
  ensembldb::select(edb, filter, columns) %>%
  as_tibble()

write.csv(tbl, "~/path/to/tbl_protein_coding_genes.csv", row.names = FALSE

########################### Vst ################
Seq1 = Americas %>% dplyr::filter(Population == "Maya")
Seq2 = Americas %>% dplyr::filter(Population == "QC_H")
test_combined_vst = rbind(Seq1, Seq2)

getVst <- function(data) {
  if ((max(data) - min(data)) > 1){
    dat1 <- data[1:(length(Seq1$Population))]
    dat2 <- data[(length(Seq1$Population)+1):(length(Seq1$Population) + (length(Seq2$Population)))]
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

full_dataset = as.data.frame(apply(test_combined_vst[3:ncol(test_combined_vst)],2, getVst))
full_dataset$Gene = row.names(full_dataset)
rownames(full_dataset) <- NULL
colnames(full_dataset) = c("VST", "Gene")

######################################### Calculate Medians #############################
Seq1_median = apply(Seq1[3:ncol(test_combined_vst)], 2, median)
Seq2_median = apply(Seq2[3:ncol(test_combined_vst)], 2, median)
median_across = as.data.frame(abs(Seq1_median - Seq2_median))
median_across$Gene = row.names(median_across)
rownames(median_across) <- NULL
colnames(median_across) = c("Medians", "Gene")

full_dataset$Medians = median_across$Medians
######################################### Filter only the Top VST Hits ##############################
# Since genes tend to have slightly different starts and ends based on how they are annotated in the UCSC genome browser 
# resulting in multiple locationsfor the same genes, we filtered to include only the genes versions that had the highest Vst value 

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
locations = read.delim("~/path/to/gene/locations.txt", header = T) # all genetic features

locations_unique = locations %>%
  dplyr::distinct(across(Gene),
                  .keep_all = TRUE)

locations_adjusted = locations_unique[locations_unique$Gene %in% filtered_dataset$Gene,]
locations_adjusted_sorted <- locations_adjusted[order(locations_adjusted$Gene),]
filtered_dataset_sorted <- filtered_dataset[order(filtered_dataset$Gene),]

############################### Filter out usually high copy numbered genes/clear errors #############
# remove  genes that have zero values as they do not have a difference in copy number and bloat the file
# Notes on why these four genes are removed are located in Supplementary methods

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

Quechua = ggplot(Outlier_removed, aes(x = Medians, y = VST, color = combined_metric)) +
  geom_point( size = 3) +
  scale_color_gradientn(
    colors = rev(viridis(256, option = "inferno")), # Reverses the color order
    oob = scales::squish) +
  ylim(0, 0.35) +
  xlim(0, 11) +
  xlab("Median Gene Copy Number Difference") + 
  ylab("VST") + 
  ggtitle("Quechua from Cerro de Pasco vs. Maya from Chiapas") + 
  labs(color = " VST x Median") +
  theme_minimal() +
  theme(legend.position = "bottom")

pdf(file = "~/path/to/Q_H_Maya_VST_Medians.pdf", width = 8, height = 9)
print(Quechua)
dev.off()
