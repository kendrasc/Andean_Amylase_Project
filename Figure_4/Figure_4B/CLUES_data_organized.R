library(RcppCNPy)
library(tidyverse)


########################## Quechua #################
prefix = "Quechua_PhasedChr1_no_missing_biallelic_no_EUR_bypop_sampled_CLUES_coal"      
epochs  = npyLoad(paste0(prefix, ".epochs.npy"))
freqs   = npyLoad(paste0(prefix, ".freqs.npy"))
logpost = npyLoad(paste0(prefix, ".post.npy"))
post    = exp(logpost)

band_frac = 1e-3
gen_max   = apply(post, 2, max)

peak_idx  = apply(post, 2, which.max)

band_idx = lapply(seq_along(gen_max), function(j) {
  which(post[, j] >= gen_max[j] * band_frac)
})

peak_freq = freqs[peak_idx]
band_low  = map_dbl(band_idx, ~ min(freqs[.x]))
band_high = map_dbl(band_idx, ~ max(freqs[.x]))

df_quechua = tibble(
  gen       = epochs[-length(epochs)],
  peak_freq = peak_freq,
  ymin      = band_low,
  ymax      = band_high
)
df_quechua$pop = "Quechua"




########################### Maya ######################
prefix = "Maya_PhasedChr1_no_missing_biallelic_no_EUR_bypop_sampled_CLUES_coal"
epochs  = npyLoad(paste0(prefix, ".epochs.npy"))
freqs   = npyLoad(paste0(prefix, ".freqs.npy"))
logpost = npyLoad(paste0(prefix, ".post.npy"))
post    = exp(logpost)

band_frac = 1e-3

gen_max   = apply(post, 2, max)

peak_idx  = apply(post, 2, which.max)

band_idx = lapply(seq_along(gen_max), function(j) {
  which(post[, j] >= gen_max[j] * band_frac)
})

peak_freq = freqs[peak_idx]
band_low  = map_dbl(band_idx, ~ min(freqs[.x]))
band_high = map_dbl(band_idx, ~ max(freqs[.x]))


df_maya = tibble(
  gen       = epochs[-length(epochs)],   # or seq_along(epochs[-1]) - 1
  peak_freq = peak_freq,
  ymin      = band_low,
  ymax      = band_high
)
df_maya$pop = "Maya"



test$pop = factor(test$pop, levels = c("Maya", "Quechua"))

cols = c("Maya" = "#800080",
          "Quechua" = "#D4AA00")

p = ggplot(test, aes(x = 28*gen)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax,
                  fill  = pop,
                  group = pop),
              alpha = 0.25, linetype = 0) +
  geom_line(aes(y = peak_freq,
                colour = pop),
            size = 1) +
  scale_fill_manual(values = cols, name = "Population") +
  scale_colour_manual(values = cols, name = "Population") +
  scale_y_continuous(limits = c(0, 0.6)) +
  geom_vline(xintercept = 6000, linetype = "dashed") +
  geom_vline(xintercept = 10000, linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(y = "SNV Frequency",
       x = "Years Before Present") +
  ggtitle("rs143597860 Frequency Over Time")

ggsave("~/Desktop/Amylase_Americas/PDFs/Quechua_Maya_rs143597860.pdf", p, width = 6, height = 4)
