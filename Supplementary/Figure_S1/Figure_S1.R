library(ggplot2)
library(tidyverse)
library(ggthemes)
library(rnaturalearth)
library(ggpubr)

# GenetoCN data
####################### GenetoCN for Amylase. #####################################################################################
genetocn = read.delim("~/Desktop/Amylase_Americas/All_amylase_gene_copy_number.txt", header=T)
OKGenomes = read.delim("~/Desktop/Amylase_Americas/1k_genomes.txt")
HGDPGenomes = read.delim("~/Desktop/Amylase_Americas/HGDP_regions.txt")

genetocn$AMY1_rounded = round(genetocn$AMY1)


EAS = c("CDX","KHV", "JPT","Tibetans", "CHB", "CHS", "Cambodian_Cambodia", "Dai_China",
        "Daur_China", "Yi_China", "Tu_China", "Tibetan", "She_China", "Mongolian_China",
        "Lahu_China", "Hezhen_China", "Yakut_Russia", "Northern_Han", "Miao_China",
        "Xibo_China", "Tujia_China", "Oroqen_China", "Naxi_China", "Japanese_Japan", "Han_China")
AFR = c("LWK", "MSL", "ACB", "ESN", "ASW", "Bantu_Kenya", "Bantu_SA", "Biaka_CAR",
        "GWD", "YRI", "Yoruba_Nigeria", "San_Namibia", "Mozabite_Algeria", "Mbuti_DRC",
        "Mandenka_Senegal")
AMR_other = c("Maya_Chiapas", "Colombian_Colombia", 
              "PUR", "MXL", "CLM", "Surui_Brazil", "Maya_Yucatan",
              "Karitiana_Brazil")
AMR_peru_pima = c("Lowland_Quechua", "Highland_Quechua", "PEL", "Pima_Mexico")
AMR = c("Maya_Chiapas", "Colombian_Colombia", "Lowland_Quechua", "Highland_Quechua",
        "PUR", "MXL", "CLM", "Surui_Brazil", "Maya_Yucatan", "PEL", "Pima_Mexico", "Karitiana_Brazil")
EUR = c("TSI", "GBR", "FIN", "CEU", "IBS", "Adygei_Russia", "Basque_France", "Bergamo_Italy",
        "French_France", "Russian_Russia", "Orcadian_Orkney", "Tuscan_Italy", "Sardinian_Italy")
SAS = c("ITU", "GIH", "PJL", "BEB", "Balochi_Pakistan", "Brahui_Pakistan", "Burusho_Pakistan",
        "Hazara_Pakistan", "Makrani_Pakistan", "Kalash_Pakistan", "Sindhi_Pakistan", "Uygur_China",
        "Pathan_Pakistan")
Middle_East = c("Bedouin_Israel", "Druze_Israel", "Palestinian_Israel", "UAE", "IRAQ", "YEMEN", "SYRIA", "SAUDI")
Oceania = c("Bougainville_Bougainville", "Papuan_PNG")

genetocn$Greater_Region = with(genetocn, 
                               ifelse(Pop %in% EAS , "EAS",
                                      ifelse(Pop %in% AFR, "AFR",
                                                    ifelse(Pop %in% EUR, "EUR", 
                                                           ifelse(Pop %in% Middle_East, "Middle East",
                                                                  ifelse(Pop %in% Oceania, "Oceania", 
                                                                         ifelse(Pop %in% AMR, "AMR", "SAS")))))))

figure_1 = genetocn %>%
  ggplot(., aes(x = Greater_Region, y = AMY1_rounded, color = Pop)) +
  geom_boxplot() +
  stat_compare_means() +
  theme_minimal() +
  theme( legend.position = "none") +
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size = 16, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

genetocn$Greater_Region_AMR_split = with(genetocn, 
                                         ifelse(Pop %in% EAS , "EAS",
                                                ifelse(Pop %in% AFR, "AFR",
                                                       ifelse(Pop %in% AMR_peru_pima, "Peruvian/Akimel O'odham",
                                                              ifelse(Pop %in% EUR, "EUR", 
                                                                     ifelse(Pop %in% Middle_East, "Middle East",
                                                                            ifelse(Pop %in% Oceania, "Oceania", 
                                                                                   ifelse(Pop %in% AMR_other, "AMR low copy Pops", "SAS"))))))))


figure_2 = genetocn %>%
  ggplot(., aes(x = Greater_Region_AMR_split, y = AMY1_rounded, color = Pop)) +
  geom_boxplot() +
  stat_compare_means() +
  theme_minimal() +
  theme( legend.position = "none") +
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size = 16, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())


combined = ggarrange(figure_1, figure_2, nrow = 2, ncol = 1, 
                     labels = c("A", "B"),
                     font.label = list(size = 24, face = "bold"),
                     label.y = 1,
                     vjust = 0.3 )


figure = annotate_figure(
  combined,
  top    = text_grob("Kruskal-Wallis tests by Geographic Region", face = "bold", size = 16),
  bottom = text_grob("Populations", size = 18),
  left   = text_grob("AMY1 Copy Number", rot = 90, size = 18)
)


pdf("~/Desktop/Amylase_Americas/PDFs/Figure_S1.pdf", width = 18, height = 8)
print(figure)
dev.off()
