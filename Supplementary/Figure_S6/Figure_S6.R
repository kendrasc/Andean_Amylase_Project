library(tidyverse)

Quechua = readRDS("~/Desktop/Amylase_Americas/LD_Quechua.hap.ld.rds", header = T)

LD_snp_of_interest_front_1 <- Quechua %>%
  dplyr::filter(POS2 == "103614521") %>%
  mutate(SNP_group = "rs143597860")

LD_snp_of_interest_back_1 <- Quechua %>%
  dplyr::filter(POS1 == "103614521") %>%
  mutate(SNP_group = "rs143597860")

LD_snp_of_interest_back_1[c(2, 3)] <- LD_snp_of_interest_back_1[c(3, 2)]
combined_derived <- rbind(LD_snp_of_interest_front_1, LD_snp_of_interest_back_1)

LD_snp_of_interest_front_2 <- Quechua %>%
  dplyr::filter(POS2 == "103825276") %>%
  mutate(SNP_group = "rs1930184")

LD_snp_of_interest_back_2 <- Quechua %>%
  dplyr::filter(POS1 == "103825276") %>%
  mutate(SNP_group = "rs1930184")

LD_snp_of_interest_back_2[c(2, 3)] <- LD_snp_of_interest_back_2[c(3, 2)]
combined_ancestral <- rbind(LD_snp_of_interest_front_2, LD_snp_of_interest_back_2)
LD_combined <- rbind(combined_derived, combined_ancestral)

figure = ggplot(LD_combined, aes(x = as.numeric(POS1), y = R.2, color = SNP_group)) +
  geom_point(alpha = 0.75) +
  ylim(0, 1) +
  scale_color_manual(  name = "Focal SNV",
                       values = c("rs143597860" = "gold", "rs1930184" = "grey")) +
  theme_classic() +
  theme(legend.position = "bottom") +
  ggtitle("R^2 for rs1930184 and rs143597860") +
  ylab("R^2") +
  xlab("Position") +
  geom_vline(xintercept = 103348464, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 103830994, linetype = "dashed", color = "black")


pdf("~/Desktop/Amylase_Americas/PDFs/Figure_S5.pdf", width = 5, height = 5)
print(figure)
dev.off()
