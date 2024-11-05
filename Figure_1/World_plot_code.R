library(ggplot2)
library(tidyverse)
library(ggthemes)
library(rnaturalearth)

# GenetoCN data
####################### GenetoCN for Amylase. #####################################################################################
genetocn = read.delim("All_amylase_gene_copy_number.txt", header=T)
OKGenomes = read.delim("1k_genomes.txt")
HGDPGenomes = read.delim("HGDP_regions.txt")

EAS = c("CDX","KHV", "JPT","Tibetans", "CHB", "CHS", "Cambodian_Cambodia", "Dai_China",
        "Daur_China", "Yi_China", "Uygur_China", "Tu_China", "Tibetan", "She_China", "Mongolian_China",
        "Lahu_China", "Hezhen_China", "Yakut_Russia", "Northern_Han", "Miao_China",
        "Xibo_China", "Tujia_China", "Oroqen_China", "Naxi_China", "Japanese_Japan", "Han_China")
AFR = c("LWK", "MSL", "ACB", "ESN", "ASW", "Bantu_Kenya", "Bantu_SA", "Biaka_CAR",
        "GWD", "YRI", "Yoruba_Nigeria", "San_Namibia", "Mozabite_Algeria", "Mbuti_DRC",
        "Mandenka_Senegal")
AMR = c("Lowland_Andeans", "Highland_Andeans", "Maya_Chiapas", "Colombian_Colombia", 
        "PUR", "MXL", "CLM", "PEL", "Pima_Mexico", "Surui_Brazil", "Maya_Mexico",
        "Karitiana_Brazil")
EUR = c("TSI", "GBR", "FIN", "CEU", "IBS", "Adygei_Russia", "Basque_France", "Bergamo_Italy",
        "French_France", "Russian_Russia", "Orcadian_Orkney", "Tuscan_Italy", "Sardinian_Italy")
SAS = c("ITU", "GIH", "PJL", "BEB", "Balochi_Pakistan", "Brahui_Pakistan", "Burusho_Pakistan",
        "Hazara_Pakistan", "Makrani_Pakistan", "Kalash_Pakistan", "Sindhi_Pakistan",
        "Pathan_Pakistan")
Middle_East = c("Bedouin_Israel", "Druze_Israel", "Palestinian_Israel", "UAE", "IRAQ", "YEMEN", "SYRIA", "SAUDI")
Oceania = c("Bougainville_Bougainville", "Papuan_PNG")

genetocn$Greater_Region = with(genetocn, 
                               ifelse(Pop %in% EAS , "EAS",
                                      ifelse(Pop %in% AFR, "AFR",
                                             ifelse(Pop %in% AMR, "AMR",
                                                    ifelse(Pop %in% EUR, "EUR", 
                                                           ifelse(Pop %in% Middle_East, "Middle_East",
                                                                  ifelse(Pop %in% Oceania, "Oceania", "SAS")))))))

uniq = unique(genetocn$Pop) # gives me a vector of all of the uniq population names

oneklocation = OKGenomes[c("Population.code", "Population.latitude", "Population.longitude")]
colnames(oneklocation) = c("Population.name", "Population.latitude", "Population.longitude")
hgdplocation = HGDPGenomes[c("Population.name", "Population.latitude", "Population.longitude")]
Population.name = c("Lowland_Andeans", "Highland_Andeans", "Maya_Chiapas", "Tibetans", "SAUDI", "YEMEN", "UAE", "SYRIA", "IRAQ")
# Saudi Arabia = Riyadh
# Yemen = Taizz
# UAE = Abu Dhabi
# Syria = Homs
# Iraq = Baghdad

Population.latitude = c(-12.04000, -10.6674, 16.7569, 30.74822, 24.633333, 13.578889, 24.466667, 34.730833, 33.315278)
Population.longitude = c(-77.030000, -76.2540, -93.1292, 92.81752, 46.716667, 44.021944, 54.366667, 36.709444, 44.366111)
additionallocation = data.frame(Population.name, Population.latitude, Population.longitude)
locations = rbind(oneklocation, hgdplocation, additionallocation)

####################### Calculate Medians
getMedians_1k <- function(data){
  pop = genetocn %>%
    dplyr::filter(Pop == data) 
  Median = median(pop$AMY1)
  return(Median)
}

test = lapply(uniq, getMedians_1k)
medians = unlist(test)
Medians_with_pop = data.frame(uniq, medians)
Medians_with_pop = Medians_with_pop[order(Medians_with_pop$uniq),c(1,2)]  %>%
  dplyr::filter(uniq != "long_read_1k" & uniq != "Colombian_Colombia" & uniq != "Bantu_SA" & uniq != "San_Namibia") # remove the small populations
OKGenomes_sorted = locations[order(locations$Population.name), c(1:3)] %>%
  dplyr::filter(Population.name != "Papuan Sepik" & Population.name != "Colombian" & Population.name != "Bantu South Africa" & Population.name != "San")# seems to work right 
OKGenomes_medians = data.frame(Medians_with_pop, OKGenomes_sorted$Population.latitude, OKGenomes_sorted$Population.longitude)
OKGenomes_medians[order(OKGenomes_medians$medians), c(1:4)]

OKGenomes_medians$Superpopulation.code = with(OKGenomes_medians, 
                                              ifelse(uniq %in% EAS , "East Asia",
                                                     ifelse(uniq %in% AFR, "Africa",
                                                            ifelse(uniq %in% AMR, "Americas",
                                                                   ifelse(uniq %in% EUR, "Europe", 
                                                                          ifelse(uniq %in% Middle_East, "Middle East",
                                                                                 ifelse(uniq %in% Oceania, "Oceania", "South Asia")))))))


world <- map_data("world")

new_world <- ggplot() +
  geom_polygon(data = world, aes(x=long, y = lat, group = group), fill="grey", alpha=0.3) +
  geom_point(data=OKGenomes_medians, aes(x=as.numeric(OKGenomes_sorted.Population.longitude), 
                                         y=as.numeric(OKGenomes_sorted.Population.latitude), 
                                         size=medians, 
                                         color=Superpopulation.code), stroke=F, alpha=0.5) +
             stroke = F, alpha = 0.5)
  scale_color_manual(values=c("South Asia" = "sienna4",
                              "Oceania" = "gray35", 
                              "Middle East" = "mediumpurple4",
                              "Europe"= "orange", 
                              "East Asia" = "maroon4", 
                              "Americas" = "darkgoldenrod3", 
                              "Africa" = "deeppink3")) +
  labs(color = "Super Population", size = "Population Median Copy Number") +
  scale_size_continuous(range=c(1,15), breaks=c( 6.1, 7, 8, 9, 10)) +
  guides(color = guide_legend(override.aes = list(size = 9)),
         fill = guide_legend(override.aes = list(size = 9))) +
  theme_void() + 
  theme(
    legend.position = "right",
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#ffffff", color = NA), 
    panel.background = element_rect(fill = "#ffffff", color = NA), 
    legend.background = element_rect(fill = "#ffffff", color = NA)
  )

new_world
pdf(file = "new_world_correct_colors_amylase_side_legend.pdf", width = 18, height = 8)
print(new_world)
dev.off()



