############################################################################
library(dplyr)
library(ggplot2)

####################### GenetoCN for Amylase. #####################################################################################
genetocn = read.delim("~/Desktop/Amylase_Americas/All_amylase_gene_copy_number.txt", header=T)
# This is the main file where I kept all of the copy  number data. The values for this are in supplementary table 1

############################################################################################################
# Geographic labels are coming from the datasets

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
EUR = c("TSI", "GBR", "FIN", "CEU", "IBS", "Adygei_Russia", "Basque_France", "Bergamo_Italy",
        "French_France", "Russian_Russia", "Orcadian_Orkney", "Tuscan_Italy", "Sardinian_Italy")
SAS = c("ITU", "GIH", "PJL", "BEB", "Balochi_Pakistan", "Brahui_Pakistan", "Burusho_Pakistan",
        "Hazara_Pakistan", "Makrani_Pakistan", "Kalash_Pakistan", "Sindhi_Pakistan", "Uygur_China",
        "Pathan_Pakistan")
Middle_East = c("Bedouin_Israel", "Druze_Israel", "Palestinian_Israel", "UAE", "IRAQ", "YEMEN", "SYRIA", "SAUDI")
Oceania = c("Bougainville_Bougainville", "Papuan_PNG")


genetocn$Greater_Region = with(genetocn, 
                               ifelse(Pop %in% EAS , "East Asia",
                                      ifelse(Pop %in% AFR, "Africa",
                                             ifelse(Pop %in% AMR_peru_pima, "PEL/Quechua/Akimel O'odham",
                                                    ifelse(Pop %in% EUR, "Europe", 
                                                           ifelse(Pop %in% Middle_East, "Middle East",
                                                                  ifelse(Pop %in% Oceania, "Oceania", 
                                                                         ifelse(Pop %in% AMR_other, "Additional Americans", "South Asia"))))))))

genetocn$AMY1_rounded = round(genetocn$AMY1)

figure = genetocn %>%
  dplyr::filter(Pop != "Colombian_Colombia" & Pop != "Bantu_SA" & Pop != "San_Namibia") %>% # These samples are filtered out for having too few individuals
  dplyr::mutate(Pop = ifelse(Pop %in% c("Highland_Quechua", "Lowland_Quechua"), "Quechua", # Quechua populations are grouped together as they should be genetically very similar
                              ifelse(Pop == "Maya_Mexico", "Maya_Yucatan", Pop))) %>%
  dplyr::mutate(Pop = dplyr::recode(Pop, "Tibetan" = "East Asia", "Cambodian_Cambodia" = "Cambodian", 
                      "Dai_China" = "Dai", "Daur_China" = "Daur", "Yi_China" = "Yi", 
                      "Uygur_China" = "Uygur", "Tu_China" = "Tu", "She_China" = "She", 
                      "Mongolian_China" = "Mongolian", "Lahu_China" = "Lahu", "Hezhen_China" = "Hezhen", 
                      "Yakut_Russia" = "Yakut", "Northern_Han" = "Northern Han", "Miao_China" = "Miao", 
                      "Xibo_China" = "Xibo", "Tujia_China" = "Tujia", "Oroqen_China" = "Oroqen", 
                      "Naxi_China" = "Naxi", "Japanese_Japan" = "Japanese", "Han_China" = "Southern Han",
                      # Africa
                      "Bantu_Kenya" = "Bantu", "Biaka_CAR" = "Biaka",
                      "Yoruba_Nigeria" = "Yoruba", "Mozabite_Algeria" = "Mozabite", "Mbuti_DRC" = "Mbuti", 
                      "Mandenka_Senegal" = "Mandenka",
                      # Americas
                      "Lowland_Andeans" = "Lowland_Andeans", "Highland_Andeans" = "Highland_Andeans", "Maya_Chiapas" = "Maya (Chiapas)", 
                      "Pima_Mexico" = "Pima", "Surui_Brazil" = "Suruí", "Maya_Yucatan" = "Maya (Yucatán)", 
                      "Karitiana_Brazil" = "Karitiana", 
                      # Europe
                      "Adygei_Russia" = "Adygei", "Basque_France" = "Basque", "Bergamo_Italy" = "Bergamo", 
                      "French_France" = "French", "Russian_Russia" = "Russian", "Orcadian_Orkney" = "Orcadian", 
                      "Tuscan_Italy" = "Tuscan", "Sardinian_Italy" = "Sardinian",
                      # South Asia
                      "Balochi_Pakistan" = "Balochi", "Brahui_Pakistan" = "Brahui", "Burusho_Pakistan" = "Burusho",
                      "Sindhi_Pakistan" = "Sindhi", "Hazara_Pakistan" = "Hazara", "Pathan_Pakistan" = "Pathan",
                      "Makrani_Pakistan" = "Makrani", "Kalash_Pakistan" = "Kalash", 
                      # Middle East
                      "Bedouin_Israel" = "Bedouin", "Druze_Israel" = "Druze", "Palestinian_Israel" = "Palestinian", 
                      "UAE" = "Emerati", "IRAQ" = "Iraqi Arab", "YEMEN" = "Yemeni", "SYRIA" = "Syrian", "SAUDI" = "Saudi",
                      # Oceania
                      "Bougainville_Bougainville" = "Bougainville", "Papuan_PNG" = "Papua New Guinean")
) %>%
  ggplot(., aes(y = reorder(Pop, -AMY1_rounded, median), x = AMY1_rounded, color=Greater_Region, fill=Greater_Region)) +
    geom_boxplot(alpha = 0.4) +
  scale_color_manual(values = c("South Asia" = "#8b4726ff",
                               "Oceania" = "#333333ff", 
                               "Middle East" = "#5d478bff",
                               "Europe" = "#ff6600ff", 
                               "East Asia" = "#ffaaeeff", 
                               "PEL/Quechua/Akimel O'odham" = "#d4aa00ff", 
                               "Additional Americans" = "#800080ff",
                               "Africa" = "#cd1076ff")) +
  scale_fill_manual(values = c("South Asia" = "#8b4726ff",
                               "Oceania" = "#333333ff", 
                               "Middle East" = "#5d478bff",
                               "Europe" = "#ff6600ff", 
                               "East Asia" = "#ffaaeeff", 
                               "PEL/Quechua/Akimel O'odham" = "#d4aa00ff", 
                               "Additional Americans" = "#800080ff",
                               "Africa" = "#cd1076ff")) +
    ylab("Population") +
    xlab("AMY1 Diploid Copy Number") +
    geom_vline(xintercept = median(genetocn %>% 
                                      dplyr::filter(Pop != "Colombian_Colombia" & Pop != "Bantu_SA" & Pop != "San_Namibia") %>% 
                                      pull(AMY1_rounded)), # add a median line
            linetype = "dashed", linewidth = 1.2) +
    guides(color=guide_legend(title="Geographic Region"), fill=guide_legend(title="Geographic Region")) +
    theme_minimal() +
    theme(legend.position = "none",
      axis.text.x = element_text(angle = 90))
 
pdf(file = "~/Desktop/Amylase_Americas/PDFs/Amylase_world_wide_boxplot.pdf", height = 14, width = 6) 
figure
dev.off()

