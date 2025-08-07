library(tidyverse)
library(ggpubr)

# --------------------------- Load & prep ---------------------------
csv_path <- "~/Desktop/Unaltered_AMY1_linked_SNPs.csv"  # change if needed
sig_snps <- read.csv(csv_path, header = TRUE)
sig_snps$Kruskal_Wallis_status_SNV_IDs <- sub("^X", "", sig_snps$Kruskal_Wallis_status_SNV_IDs)

# Population order (includes all you referenced; extra pops appended if present)
pop_order_base <- c(
  "PUR","MXL","CLM","Maya","PEL","Quechua",      # American (Quechua–PUR)
  "JPT","KHV","CHB","CHS","CDX",                 # East Asia (CDX–JPT)
  "PJL","ITU","GIH","BEB","STU",                 # South Asia (PJL–ITU)
  "GBR","FIN","CEU","IBS","TSI",                 # Europe (GBR–FIN)
  "ASW","YRI","LWK","ESN","GWD","MSL","ACB"      # Africa (ASW–YRI)
)
pop_order <- unique(c(pop_order_base, setdiff(unique(sig_snps$Pop), pop_order_base)))

# ------------------- Alpha (piecewise; robust) --------------------
# < 0.01 => strong alpha (0.75–1.00) ; 0.01–0.05 => faint (0.05–0.30)
compute_alpha <- function(p) {
  p_clamped <- pmin(pmax(p, 1e-100), 0.05)
  nl10 <- -log10(p_clamped)
  nl10_lo  <- 1.30103   # -log10(0.05)
  nl10_mid <- 2.00000   # -log10(0.01)
 
  nl10_max_lo <- suppressWarnings(max(nl10[p_clamped < 0.01], na.rm = TRUE))
  if (!is.finite(nl10_max_lo)) nl10_max_lo <- nl10_mid
  if (nl10_max_lo <= nl10_mid) nl10_max_lo <- nl10_mid + 1e-6
 
  denom_hi <- (nl10_mid - nl10_lo)         # 0.69897
  denom_lo <- (nl10_max_lo - nl10_mid)
 
  alpha_val <- ifelse(
    p_clamped < 0.01,
    0.75 + 0.25 * (nl10 - nl10_mid) / denom_lo,            # -> [0.75, 1.00]
    0.05 + (nl10 - nl10_lo) / denom_hi * (0.30 - 0.05)     # -> [0.05, 0.30]
  )
  pmin(pmax(alpha_val, 0.01), 1)
}

# ------------------------- Data frames ----------------------------
# All significant points
base_df <- sig_snps %>%
  filter(!is.na(KW_adj_values), KW_adj_values < 0.05) %>%
  transmute(
    SNV_ID  = Kruskal_Wallis_status_SNV_IDs,
    SNV_pos = as.numeric(Kruskal_Wallis_status_SNV_IDs),
    Pop     = factor(Pop, levels = pop_order),
    KW_adj_values,
    alpha_val = compute_alpha(KW_adj_values)
  ) %>%
  arrange(alpha_val) %>%
  mutate(alpha_base = pmax(alpha_val * 0.6, 0.03))  # damp base so highlights pop

# Highlights: Quechua & PEL (stronger alpha floor)
hl_df <- base_df %>%
  filter(Pop %in% c("Quechua","PEL")) %>%
  mutate(
    alpha_hl = ifelse(KW_adj_values < 0.01,
                      pmax(alpha_val, 0.85),
                      pmax(alpha_val, 0.35))
  )

# --------------------- Continent separators (lines only) ----------
regions_sets <- list(
  American     = c("Quechua","PEL","Maya","CLM","MXL","PUR"),
  "East Asia"  = c("CDX","KHV","JPT","CHB","CHS"),
  "South Asia" = c("PJL","ITU","GIH","BEB","STU"),
  Europe       = c("GBR","FIN","CEU","IBS","TSI"),
  Africa       = c("ASW","YRI","LWK","ESN","GWD","MSL","ACB")
)
ypos <- setNames(seq_along(pop_order), pop_order)
region_spans <- purrr::map_dfr(names(regions_sets), function(reg){
  pops <- intersect(pop_order, regions_sets[[reg]])
  if (length(pops) == 0) return(NULL)
  idx <- sort(unname(ypos[pops]))
  tibble(region = reg, ymin = min(idx) - 0.5, ymax = max(idx) + 0.5)
}) %>% arrange(ymin)
div_lines <- if (nrow(region_spans) > 1) region_spans$ymax[-nrow(region_spans)] else numeric(0)

# --------------------------- Plot --------------------------------
xmin <- 103347535; xmax <- 103937847
highlight_cols <- c(Quechua = "#D81B60", PEL = "#1E88E5")

p <- ggplot() +
  # Thin divider lines between regions
  (if (length(div_lines) > 0)
    geom_hline(yintercept = div_lines, color = "grey60", linewidth = 0.4)
   else NULL) +
 
  xlim(xmin, xmax) +
  geom_vline(xintercept = c(103348464, 103830994), linetype = "dashed") +
 
  # All significant SNVs in grey, faded by alpha
  geom_tile(
    data = base_df,
    aes(x = SNV_pos, y = Pop, alpha = alpha_base),
    width = 1500, height = 0.8,
    fill = "grey60", colour = NA
  ) +
 
  # Emphasize Quechua and PEL (color + alpha floor)
  geom_tile(
    data = hl_df %>% filter(Pop == "Quechua"),
    aes(x = SNV_pos, y = Pop, alpha = alpha_hl),
    width = 1500, height = 0.8,
    fill = highlight_cols["Quechua"], colour = NA
  ) +
  geom_tile(
    data = hl_df %>% filter(Pop == "PEL"),
    aes(x = SNV_pos, y = Pop, alpha = alpha_hl),
    width = 1500, height = 0.8,
    fill = highlight_cols["PEL"], colour = NA
  ) +
 
  labs(
    x = "Genomic position (bp)",
    y = "Population",
    title = "SNVs Linked to AMY1 Copy Number (Bonferroni < 0.05)",
    subtitle = "Continents separated by thin lines; Quechua & PEL emphasized"
  ) +
  scale_y_discrete(drop = FALSE) +
  scale_alpha_identity() +
  theme_classic() +
  theme(legend.position = "none")

p

# install once if needed
if (!requireNamespace("ragg", quietly = TRUE)) install.packages("ragg")

ggsave(
  filename = "SNVs_AMY1_continent_lines.png",
  plot     = p,
  width    = 8,       # inches
  height   = 5,        # increase if you want more vertical room for populations
  units    = "in",
  dpi      = 600,      # print quality
  device   = ragg::agg_png,
  bg       = "white"
)


#############LEGEND

# Standalone legend for alpha mapping and save as PNG
# (Uses cowplot to extract the legend and ragg for crisp output)

if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("cowplot", quietly = TRUE)) install.packages("cowplot")
if (!requireNamespace("ragg", quietly = TRUE)) install.packages("ragg")

library(ggplot2)
library(cowplot)

# Dummy data just to generate the legend
legend_df <- data.frame(
  band = factor(c("< 0.01", "0.01–0.05"), levels = c("< 0.01", "0.01–0.05")),
  x = 1, y = 1
)

# Build a plot that only serves to create the legend
p_for_legend <- ggplot(legend_df, aes(x = x, y = y, alpha = band)) +
  geom_point(shape = 15, size = 8, color = "grey30") +
  scale_alpha_manual(
    name   = "Adjusted p (Bonf.)",
    values = c("< 0.01" = 0.95, "0.01–0.05" = 0.20)
  ) +
  guides(alpha = guide_legend(override.aes = list(shape = 15, size = 8, color = "grey30"))) +
  theme_void() +
  theme(legend.position = "right")

# Extract legend as a grob and wrap in a draw plot
leg <- cowplot::get_legend(p_for_legend)
legend_plot <- cowplot::ggdraw(leg)

# Save high-resolution PNG
ggsave(
  filename = "alpha_legend.png",
  plot     = legend_plot,
  width    = 2.4,   # inches; adjust to taste
  height   = 1.4,   # inches
  units    = "in",
  dpi      = 600,
  device   = ragg::agg_png,
  bg       = "white"
)
