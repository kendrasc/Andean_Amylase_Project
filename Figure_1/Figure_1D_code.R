############################################################################
library(dplyr)
library(ggplot2)
library(forcats)

####################### GenetoCN for Amylase. #####################################################################################
genetocn = read.delim("~/Desktop/Amylase_Americas/All_amylase_gene_copy_number.txt", header=T)

############################################################################################################
AMR = c("Lowland_Quechua", "Highland_Quechua", "Maya_Chiapas",
        "PUR", "MXL", "CLM", "PEL", "Pima_Mexico", "Surui_Brazil", "Maya_Yucatan",
        "Karitiana_Brazil")
high = c("Pima_Mexico", "Lowland_Quechua", "Highland_Quechua", "PEL")
genetocn$Copy_Number = with(genetocn,
                            ifelse (Pop %in% high, "Yes", "No"))

########################### Create Plot ###########################
Americas = genetocn %>%
  dplyr::filter(Pop %in% AMR) %>%
  dplyr::mutate(
    Pop = dplyr::recode(Pop, "Highland_Quechua" = "Quechua",
                        "Lowland_Quechua" = "Quechua")) %>%
  dplyr::mutate(Pop= fct_relevel(Pop, "Surui_Brazil", "Karitiana_Brazil","Maya_Yucatan",
                                 "Maya_Chiapas", "MXL", "PUR", "CLM",
                                 "Pima_Mexico", "Quechua", "PEL",)) %>%
  ggplot(., aes(x = round(AMY1), fill=Copy_Number, color=Copy_Number)) +
  facet_grid(rows = vars(Pop)) +
geom_vline(xintercept = seq(1, 21, by = 2), linetype = "dashed") +
  geom_histogram(alpha = 0.5, binwidth = 1) + 
  xlab("") +
  ylab("") +
  scale_color_manual(values = c("Yes" = "#d4aa00ff", "No" = "#800080ff")) +
  scale_fill_manual( values = c("Yes" = "#d4aa00ff", "No" = "#800080ff")) +
  theme_bw() +
  labs(color = "AMY1 Increase", fill = "AMY1 Increase") +
  theme( 
    legend.position="bottom",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.background = element_rect(fill = "white", color="white"),
    strip.background = element_rect(fill = "white", color = "white")) + 
  ggtitle("") 

pdf("~/Desktop/Amylase_Americas/PDFs/genetocn_amy1_AMR_rounded_histogram.pdf", height = 10, width = 8)
Americas
dev.off()

