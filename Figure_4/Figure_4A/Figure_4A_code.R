# Code to create XP-EHH 
library(ggplot2)
library(dplyr)

selscan_xpehh_test_norm = read.delim("~/Desktop/Amylase_Americas/Quechua_Maya_filtered_2000000_wagh_correct_locations.xpehh.out.norm")

selscan_xpehh_test_norm$normxpehh_no_na = selscan_xpehh_test_norm$normxpehh
selscan_xpehh_test_norm$normxpehh_no_na[is.na(selscan_xpehh_test_norm$normxpehh_no_na)] = 0

threshold_selscan_norm = quantile(abs(selscan_xpehh_test_norm$normxpehh_no_na), 0.99)
selscan_xpehh_test_norm$in_99th_percentile = ifelse(abs(selscan_xpehh_test_norm$normxpehh_no_na) >= threshold_selscan_norm, "yes", "no")

selscan_norm = selscan_xpehh_test_norm %>%
  dplyr::filter(pos > 103000000 & pos < 104300000) %>%
  ggplot(aes(pos, normxpehh_no_na, color = in_99th_percentile, alpha = 0.5)) + 
  theme_minimal() +
  geom_point() + 
  geom_hline(yintercept = quantile(abs(selscan_xpehh_test_norm$normxpehh_no_na), 0.99), 
             linetype = "dashed") +
  geom_hline(yintercept = -quantile(abs(selscan_xpehh_test_norm$normxpehh_no_na), 0.99), 
             linetype = "dashed") +
  scale_color_manual(values = c("yes" = "#a02c5aff", "no" = "#ffaaeeff")) +
  geom_rect(aes(xmin = 103626688, xmax = 103759084, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE,
            fill = "#cd1076ff" , alpha = 0.01) +
  xlab("Position (BP)") +
  ylab("XP-EHH") +
  ylim(-7,7) +
  theme(legend.position = "top") +
  ggtitle("XP-EHH in Selscan Normalized")

pdf("~/Desktop/Amylase_Americas/selscan_norm_wagh_correct_positions.pdf", width = 6, height = 5)
print(selscan_norm)
dev.off()
