############################################################################
library(dplyr)
library(ggplot2)

####################### GenetoCN for Amylase. #####################################################################################
genetocn = read.delim("~/Desktop/Amylase_Americas/All_amylase_gene_copy_number.txt", header=T)

############################################################################################################
EAS = c("CDX","KHV", "JPT","Tibetan", "CHB", "CHS", "Cambodian_Cambodia", "Dai_China",
        "Daur_China", "Yi_China", "Uygur_China", "Tu_China", "Tibetan", "She_China", "Mongolian_China",
        "Lahu_China", "Hezhen_China", "Yakut_Russia", "Northern_Han", "Miao_China",
        "Xibo_China", "Tujia_China", "Oroqen_China", "Naxi_China", "Japanese_Japan", "Han_China",
        "Dai_China")
AFR = c("LWK", "MSL", "ACB", "ESN", "ASW", "Bantu_Kenya", "Bantu_SA", "Biaka_CAR",
        "GWD", "YRI", "Yoruba_Nigeria", "San_Namibia", "Mozabite_Algeria", "Mbuti_DRC",
        "Mandenka_Senegal")
AMR = c("Lowland_Andeans", "Highland_Andeans", "Maya_Chiapas", "Colombian_Colombia", 
        "PUR", "MXL", "CLM", "PEL", "Pima_Mexico", "Surui_Brazil", "Maya_Yucatan",
        "Karitiana_Brazil", "Quechua")
EUR = c("TSI", "GBR", "FIN", "CEU", "IBS", "Adygei_Russia", "Basque_France", "Bergamo_Italy",
        "French_France", "Russian_Russia", "Orcadian_Orkney", "Tuscan_Italy", "Sardinian_Italy")
SAS = c("ITU", "GIH", "PJL", "BEB", "Balochi_Pakistan", "Brahui_Pakistan", "Burusho_Pakistan",
        "Hazara_Pakistan", "Makrani_Pakistan", "Kalash_Pakistan", "Sindhi_Pakistan",
        "Pathan_Pakistan")
Middle_East = c("Bedouin_Israel", "Druze_Israel", "Palestinian_Israel", "UAE", "IRAQ", "YEMEN", "SYRIA", "SAUDI")
Oceania = c("Bougainville_Bougainville", "Papuan_PNG")

genetocn$Greater_Region = with(genetocn, 
                               ifelse(Pop %in% EAS , "East Asia",
                                      ifelse(Pop %in% AFR, "Africa",
                                             ifelse(Pop %in% AMR, "Americas",
                                                    ifelse(Pop %in% EUR, "Europe", 
                                                           ifelse(Pop %in% Middle_East, "Middle East",
                                                                  ifelse(Pop %in% Oceania, "Oceania", "South Asia")))))))
genetocn$American = with(genetocn,
                         ifelse (Pop %in% AMR, "Indigenous American", "Other Population"))

# For AMY2A and AMY2B
figure = genetocn %>%
  dplyr::filter(Pop != "Colombian_Colombia" & Pop != "Bantu_SA" & Pop != "San_Namibia") %>%
  dplyr::mutate(Pop = ifelse(Pop %in% c("Highland_Quechua", "Lowland_Quechua"), "Quechua",
                             ifelse(Pop == "Maya_Mexico", "Maya_Yucatan", Pop))) %>%
  dplyr::mutate(Pop = dplyr::recode(Pop, "Tibetan" = "Tibetan", "Cambodian_Cambodia" = "Cambodian", 
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
  ggplot(., aes(y = reorder(Pop, -round(AMY2B), median), x = round(AMY2B), color=American, fill=American)) +
  geom_boxplot(alpha = 0.4) +
  scale_color_manual(values = c("Indigenous American" = "#d4aa00ff", "Other Population" = "#800080ff")) +
  scale_fill_manual(values = c("Indigenous American" = "#d4aa00ff", "Other Population" = "#800080ff")) +
  ylab("Population") +
  xlab("AMY2B Diploid Copy Number") +
  geom_vline(xintercept = mean(round(genetocn %>% 
                                       dplyr::filter(Pop != "Colombian_Colombia" & Pop != "Bantu_SA" & Pop != "San_Namibia") %>% 
                                       pull(AMY2B))), 
             linetype = "dashed", linewidth = 1.2) +
  guides(color=guide_legend(title="Geographic Region"), fill=guide_legend(title="Geographic Region")) +
  theme_minimal() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 90))

pdf("~/Desktop/Amylase_Americas/PDFs/World_AMY2B.pdf", width = 7, height = 10)
print(figure)
dev.off()

