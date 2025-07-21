library(tidyr)
library(ggplot2)
library(forcats)

Admixture = read.csv("~/Desktop/Amylase_Americas/CSVs/Admixture.csv")

# Reshape data to long format
Admixture_long <- Admixture %>%
  dplyr::filter(Pop != "QC_L") %>%
  pivot_longer(
    cols = starts_with("Fraction"),
    names_to = "Ancestry",
    values_to = "Value"
  )

Admixture_long$Pop <- fct_relevel(Admixture_long$Pop, "QC_H", "MAYA", "PEL", "MXL")

# Create the faceted bar plot
test = ggplot(Admixture_long, aes(x = factor(ID), y = Value, fill = Ancestry)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Pop, scales = "free_x", nrow = 1) +
  labs(
    x = "Individuals",
    y = "Fraction",
    title = "Ancestry Proportions",
    fill = "Ancestry"
  ) +
  scale_fill_manual(
    name = "Ancestry",
    labels = c("African", "Amerindian", "East Asian", "European"),
    values = c("Fraction.AFR" = "#cd1076ff", "Fraction.AMR" = "#ffcc00ff", 
               "Fraction.EAS" = "#ffaaeeff", "Fraction.EUR" = "#ff6600ff")
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(), # Remove x-axis text labels
    axis.ticks.x = element_blank(), # Remove x-axis ticks
    panel.spacing = unit(0.5, "lines") # Adjust space between facets
  )

test

pdf("~/Desktop/Amylase_Americas/PDFs/test_admixture_pdf.pdf", width = 15, height = 3)
print(test)
dev.off()


