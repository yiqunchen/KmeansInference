# ---------------------------------------------------------------------------
# plot_datathin_sigma_box.R -- SUPPLEMENTARY. Replaces the variance block of the
# data-thinning table with a boxplot: the split-variance estimates inflate from
# the true sigma=1 toward ~3 as the clusters separate, which is what "poisons"
# the data-thinning split. (Detection and conditional power are already curves in
# datathin.png; only this variance-inflation block was table-only.)
#
# Grammar: colour = variance handling (sigma_MED green, sigma_sample vermillion);
# sigma=1 (truth) = grey dashed reference line.
#
# Usage: Rscript sims/plot_datathin_sigma_box.R
# ---------------------------------------------------------------------------
source("sims/house_style.R")
suppressMessages({ library(dplyr); library(tidyr) })

d <- as.data.frame(readRDS("sims/results/datathin_compare.rds"))

long <- d %>%
  select(delta, sig_med, sig_samp) %>%
  pivot_longer(c(sig_med, sig_samp), names_to = "est", values_to = "sigma") %>%
  mutate(est   = factor(sub("^sig_", "", est), levels = c("med", "samp")),
         dfac  = factor(delta, levels = sort(unique(delta))))

fig <- ggplot(long, aes(x = dfac, y = sigma, fill = est)) +
  geom_hline(yintercept = 1, linewidth = 0.5, linetype = "dashed", colour = km_ref) +
  geom_boxplot(position = position_dodge(width = 0.78), width = 0.66,
               outlier.size = 0.35, outlier.alpha = 0.35, linewidth = 0.35,
               colour = "black") +
  scale_fill_manual(values = c(med = km_col[["med"]], samp = km_col[["sample"]]),
                    labels = c(med = "sigma-MED split", samp = "sigma-sample split"),
                    name = NULL) +
  labs(x = expression(separation~~delta),
       y = expression(hat(sigma)~~"used to split")) +
  theme_km()

ggsave_km(fig, "sims/results/datathin_sigma_box", width = 7.2, height = 4.4)
cat("Wrote datathin_sigma_box.{pdf,png}\n")
