library(tidyverse)
library(ape)
library(ggtree)
library(ggpubr)

################## PBS
# PBS calculation is based on https://mccoy-lab.github.io/hgv_modules/index.html
fst_results <- read.table("~/Desktop/Amylase_Americas/fst_results.txt", header = TRUE)

fst_results$fst.mxl.pel_no_na = fst_results$fst.mxl.pel
fst_results$fst.mxl.pel_no_na[is.na(fst_results$fst.mxl.pel_no_na) |
                                fst_results$fst.mxl.pel_no_na < 0] <- 0


fst_results$fst.mxl.yri_no_na = fst_results$fst.mxl.yri
fst_results$fst.mxl.yri_no_na[is.na(fst_results$fst.mxl.yri_no_na) |
                                fst_results$fst.mxl.yri_no_na < 0] <- 0

fst_results$fst.pel.yri_no_na = fst_results$fst.pel.yri
fst_results$fst.pel.yri_no_na[is.na(fst_results$fst.pel.yri_no_na) |
                                fst_results$fst.pel.yri_no_na < 0] <- 0

pbs <- fst_results %>%
  # calculate branch lengths between populations
  mutate(T_mxl_pel = -log(1 - fst.mxl.pel_no_na),
         T_mxl_yri = -log(1 - fst.mxl.yri_no_na),
         T_pel_yri = -log(1 - fst.pel.yri_no_na)) %>%
  # calculate pbs
  mutate(pbs = ((T_mxl_pel + T_pel_yri) - (T_mxl_yri)) / 2) %>%
  arrange(-pbs)

# Places where the value is negative just means where MXL and PEL don't really differ
pbs$pbs_no_na = pbs$pbs
pbs$pbs_no_na[is.na(pbs$pbs_no_na) |
                pbs$pbs_no_na < 0] <- 0


threshold_norm <- quantile(abs(pbs$pbs_no_na), 0.99)
pbs$in_99th_percentile <- ifelse(abs(pbs$pbs_no_na) >= threshold_norm, "yes", "no")


amy_pbs = pbs %>%
  dplyr::filter(POS == 103614521) %>%
  pull(pbs_no_na)


pbs_comparison = ggplot(data = pbs,
       aes(x = pbs_no_na, fill = in_99th_percentile)) +
 geom_histogram() +
  geom_vline(xintercept= amy_pbs, linetype = "dashed", linewidth = 2.5, color = "#a02c5aff") +
  scale_fill_manual(values = c("yes" = "#a02c5aff", "no" = "#ffaaeeff")) +
  ylab("Count") +
  xlab("PBS") +
  ggtitle("PBS Chr 1 SNPs: PEL, MXL, & YRI") +
  theme_minimal() + theme(legend.position = "bottom")



###################### XP-EHH

############## Selscan Normalized
selscan_xpehh_test_norm = read.delim("/Users/panda_bear/Desktop/Amylase_Americas/IHS/Quechua_Maya_filtered_2000000_wagh_correct_locations.xpehh.out.norm")

# Replace NA values with 0
selscan_xpehh_test_norm$normxpehh_no_na <- selscan_xpehh_test_norm$normxpehh
selscan_xpehh_test_norm$normxpehh_no_na[is.na(selscan_xpehh_test_norm$normxpehh_no_na)] <- 0

# Compute the 99th percentile threshold
threshold_selscan_norm <- quantile(abs(selscan_xpehh_test_norm$normxpehh_no_na), 0.99)
selscan_xpehh_test_norm$in_99th_percentile <- ifelse(abs(selscan_xpehh_test_norm$normxpehh_no_na) >= threshold_selscan_norm, "yes", "no")

amy_xpehh = selscan_xpehh_test_norm %>%
  dplyr::filter(pos == 103614521) %>%
  pull(normxpehh_no_na)


xpehh_comparison = ggplot(data = selscan_xpehh_test_norm,
                        aes(x = abs(normxpehh_no_na), fill = in_99th_percentile)) +
  geom_histogram() +
  geom_vline(xintercept= amy_xpehh, linetype = "dashed", linewidth = 2.5, color = "#a02c5aff") +
  scale_fill_manual(values = c("yes" = "#a02c5aff", "no" = "#ffaaeeff")) +
  ylab("Count") +
  xlab("XP-EHH") +
  ggtitle("XP-EHH Chr 1 SNPs: Quechua & Maya") +
  theme_minimal() + theme(legend.position = "bottom")



########################## Ohana calculations

map <- read.table("~/Desktop/merged_chr1.map", header = FALSE, stringsAsFactors = FALSE)
colnames(map) <- c("chr", "snp_id", "genetic_dist", "pos")

scan = read.delim("~/Desktop/lle-ratios.txt")

if (nrow(scan) != nrow(map)) {
  stop("Error: number of SNPs in map and selection scan don't match!")
}
scan_annotated <- bind_cols(map, scan)

# Replace NA values with 0
scan_annotated$lle.ratio_no_na <- scan_annotated$lle.ratio
scan_annotated$lle.ratio_no_na[is.na(scan_annotated$lle.ratio_no_na)] <- 0

# Compute the 99th percentile threshold
threshold_ohana <- quantile(scan_annotated$lle.ratio_no_na, 0.99)
scan_annotated$in_99th_percentile <- ifelse(abs(scan_annotated$lle.ratio_no_na) >= threshold_ohana, "yes", "no")

amy_ohana = scan_annotated %>%
  dplyr::filter(pos == 103614521) %>%
  pull(lle.ratio_no_na)

ohana_comparison = ggplot(data = scan_annotated,
                          aes(x = abs(lle.ratio_no_na), fill = in_99th_percentile)) +
  geom_histogram() +
  geom_vline(xintercept= amy_ohana, linetype = "dashed", linewidth = 2.5, color = "#a02c5aff") +
  scale_fill_manual(values = c("yes" = "#a02c5aff", "no" = "#ffaaeeff")) +
  ylab("Count") +
  xlab("Ohana") +
  ggtitle("Ohana Chr 1 SNVs: Quechua") +
  theme_minimal() + theme(legend.position = "bottom")


SNV_plot = ggarrange(pbs_comparison, xpehh_comparison, ohana_comparison, 
          labels = c("A", "B", "C"),
          nrow = 3, ncol = 1, common.legend = TRUE)

ggsave("~/Desktop/Amylase_Americas/PDFs/Figure_S13_new.pdf", SNV_plot, width = 6, height = 9)

