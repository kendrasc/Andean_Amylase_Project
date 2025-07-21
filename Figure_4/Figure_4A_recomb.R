##### Code to  make recombination graph
library(ggplot2)
library(tidyr)

recomb_test = read.delim("~/Desktop/chr1_avg.bedgraph", header = F)

figure = recomb_test %>%
  dplyr::filter(V1 == "chr1") %>%
  ggplot(., aes(x=V2, y = V4)) +
  geom_line() +
  xlim(103000000, 104300000) +
  ylab("Recombination Rate") +
  xlab("Position (BP)") +
  ylim(0, 225) +
  theme_minimal()

figure

pdf("~/Desktop/Amylase_Americas/PDFs/recomb_rates.pdf", width = 6, height = 2.5)
print(figure)
dev.off()

