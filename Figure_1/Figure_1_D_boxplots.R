############################################################################
library(dplyr)
library(ggplot2)

####################### GenetoCN for Amylase. #####################################################################################
genetocn = read.delim("All_amylase_gene_copy_number.txt", header=T)

###################### Relabel American Populations #################
AMR = c("Lowland_Andeans", "Highland_Andeans", "Maya_Chiapas", "Colombian_Colombia", 
        "PUR", "MXL", "CLM", "PEL", "Pima_Mexico", "Surui_Brazil", "Maya_Mexico",
        "Karitiana_Brazil")

genetocn$Copy_Number = with(genetocn,
                            ifelse (Pop %in% high, "Yes", "No"))
genetocn$American = with(genetocn,
                            ifelse (Pop %in% AMR, "Indigenous American", "Other Population"))

# Create figure and filter out populations that have low amounts of individuals
figure = genetocn %>%
  dplyr::filter(Pop != "Colombian_Colombia" & Pop != "Bantu_SA" & Pop != "San_Namibia") %>%
  dplyr::mutate(Pop = ifelse(Pop %in% c("Highland_Andeans", "Lowland_Andeans"), "Quechua",
                              ifelse(Pop == "Maya_Mexico", "Maya_Yucatan", Pop))) %>%
  ggplot(., aes(y = reorder(Pop, -round(AMY1), median), x = round(AMY1), color=American, fill=American)) +
    geom_boxplot(alpha = 0.4) +
    scale_color_manual(values = c("Indigenous American" = "darkgoldenrod3", "Other Population" = "maroon4")) +
    scale_fill_manual(values = c("Indigenous American" = "darkgoldenrod3", "Other Population" = "maroon4")) +
    xlab("Population") +
    ylab("AMY1 Copy Number") +
    geom_vline(xintercept = mean(round(genetocn %>% 
                                      dplyr::filter(Pop != "Colombian_Colombia" & Pop != "Bantu_SA" & Pop != "San_Namibia") %>% 
                                      pull(AMY1))), 
            linetype = "dashed", linewidth = 1.2) +
    guides(color=guide_legend(title="Population of Interest"), fill=guide_legend(title="Population of Interest")) +
    theme_classic() +
    theme(legend.position = "top",
      axis.text.x = element_text(angle = 90))
 
pdf(file = "Amylase_world_wide_boxplot.pdf", height = 20, width = 6) 
figure
dev.off()

