setwd("~/Desktop/Gene_CNVs")
library(dplyr)
library(ggplot2)
library(ggrepel)
library(AnnotationHub)
library(ensembldb)
library(dplyr)

######################################### Filtering data for Protein Coding Genes #######################################################################


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

######### Pull Only Protein Coding Genes (23069) #######################
supportedFilters()
filter = ~ gene_name %in% keys & gene_biotype == "protein_coding"
tbl =
  ensembldb::select(edb, filter, columns) %>%
  as_tibble()

##### Protein Coding Genes shared between the data sets (18928) ########
filtered_data <- modern %>%
  dplyr::select(any_of(tbl$GENENAME))

##### Merging the IDs with the filtered data set #######################
filtered_with_samples = data.frame(modern$X, modern$continent,
                                   modern$pop, filtered_data)


######################################### Calculating VST ################################################################################################
Seq1 = filtered_with_samples %>% dplyr::filter(modern.pop == "PEL")
Seq2 = filtered_with_samples %>% dplyr::filter(modern.pop == "CLM")
test_combined_vst = rbind(Seq1, Seq2)

getVst <- function(data) {
  if ((max(data) - min(data)) > 1){
    dat1 <- data[1:(length(Seq1$modern.pop))]
    dat2 <- data[(length(Seq1$modern.pop)+1):(length(Seq1$modern.pop) + (length(Seq2$modern.pop)))]
    Vtotal <- var(c(dat1,dat2))
    Vgroup <- ((var(dat1)*length(dat1)) + (var(dat2)*length(dat2))) /
      (length(dat1)+length(dat2))
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

full_dataset = as.data.frame(apply(test_combined_vst[4:18931],2, getVst))
full_dataset$Gene = row.names(full_dataset)
rownames(full_dataset) <- NULL
colnames(full_dataset) = c("VST", "Gene")

######################################### VST Figure #####################################################################################################
locations = read.delim("locations_please_work.txt", header = T) # all genetic features (rename)

# file to get unique VST values
full_dataset_unique = full_dataset %>%
  dplyr::distinct(across(Gene),
                  .keep_all = TRUE)

# file to get unique locations --> Choosing the first location
locations_unique = locations %>%
  dplyr::distinct(across(Gene),
                  .keep_all = TRUE)

locations_adjusted = locations_unique[locations_unique$Gene %in% full_dataset_unique$Gene,]

no_dates_dataset = full_dataset_unique %>%
  dplyr::filter(!(Gene %in% dates))

no_dates_dataset_sorted <- no_dates_dataset[order(no_dates_dataset$Gene),]
locations_adjusted_sorted <- locations_adjusted[order(locations_adjusted$Gene),]

genesOfInterest = c( "AMY1C", "AMY1A", "AMY1B")
str(genesOfInterest)

tophits = figure %>% dplyr::filter(VST > .2)
tophits_sorted = tophits[order(-tophits$VST), c(1:5)]
############################################ Create Manhattan plot ######################################

don <- figure %>% 
  # Compute chromosome size
  group_by(CHROM) %>% 
  summarise(chr_len=max(START)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len)  %>%
  dplyr::select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(figure, ., by=c("CHROM"="CHROM")) %>%
  # Add a cumulative position of each SNP
  arrange(CHROM, START) %>%
  mutate( STARTcum=START+tot) %>%
  # Add highlight and annotation information
  mutate( is_highlight=ifelse(Gene %in% genesOfInterest, "yes", "no"))

don_no_m = don %>%
  dplyr::filter(CHROM != 23)

# Prepare X axis
axisdf <- don_no_m %>% group_by(CHROM) %>% summarize(center=( max(STARTcum) + min(STARTcum) ) / 2 )

# Fix labels for X and Y
axisdf$CHROM <- as.character(axisdf$CHROM)
axisdf$CHROM[axisdf$CHROM == "24"] <- "X"
axisdf$CHROM[axisdf$CHROM == "25"] <- "Y"

# Rename one of the AMY1s (in this case, AMY1C) as "AMY1"
don_no_m$Gene[don_no_m$Gene == "AMY1C"] <- "AMY1"

# Make the plot
p <- ggplot(don_no_m, aes(x=STARTcum, y=VST, label=Gene)) +
  geom_point( aes(color=as.factor(CHROM)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("pink", "maroon"), 22 )) +
  geom_hline(yintercept = 0.2, linetype = "dashed") +
  geom_label_repel(data= subset(don_no_m, Gene %in% genesOfInterest & VST > .2), aes(label=Gene),
                   max.overlaps = Inf) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$CHROM, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) + 
  ylim(0,.3) + xlab("Chromosome") + ggtitle("Peruvians in Lima vs. Colombians in Medellín") +
  
  # Add highlighted points
  geom_point(data=subset(don, is_highlight=="yes"), color="gold", size=3) +
  geom_point(data=subset(don, is_highlight=="yes"),shape = 1,size = 3,colour = "black") +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    text  = element_text(size= 10, color = "darkred"),
    panel.background = element_rect(fill = "white", color="white"),  legend.key = element_rect(fill = "white"), plot.title=element_text( size = 15), axis.title.y = element_text( size = 12),
    axis.title.x = element_text( size = 12), axis.text.y = element_text(colour = "black"), axis.text.x = element_text(color = "black"))

pdf(file = "PEL_CLM.pdf", width = 12, height =4)
print(p)
dev.off()
