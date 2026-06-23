# ---------------------------------------------------------------------------
# plot_cw_power_box.R -- Table 2 (known-sigma power: union vs path) AS A FIGURE.
# Boxplots of the per-replicate selective p-value, by separation delta, faceted
# by q. The fraction of each box BELOW the 0.05 line is the (selective) power;
# the union boxes sit lower than the path boxes -> uniformly more powerful.
#
# Grammar: known-sigma is monochrome (oracle/black); path/union is the
# open-vs-filled distinction (path = hollow, union = filled grey), mirroring the
# point-geom shape exception in house_style.R. 0.05 = reference line (grey dashed).
#
# Usage: Rscript sims/plot_cw_power_box.R [results_dir]
#   default dir = sims/results/sweep_cw_uv  (re-point to sweep_cw_uv_exact after
#   the exact backfill completes).
# ---------------------------------------------------------------------------
source("sims/house_style.R")
suppressMessages({ library(dplyr); library(tidyr) })

DIR <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "sims/results/sweep_cw_uv"

# ---- load every per-replicate power cell ----------------------------------
cells <- list.files(file.path(DIR, "cells"), pattern = "^cw_power_q[0-9]+_d[0-9]+__c[0-9]+\\.rds$",
                    full.names = TRUE)
stopifnot(length(cells) > 0)
raw <- do.call(rbind, lapply(cells, function(f) {
  d <- readRDS(f)
  keep <- c("p_path", "p_union", "label")
  d[, keep, drop = FALSE]
}))
# parse q and delta from the label (cw_power_q<Q>_d<D>)
m <- regmatches(raw$label, regexec("cw_power_q([0-9]+)_d([0-9]+)", raw$label))
raw$q     <- as.integer(vapply(m, `[`, "", 2))
raw$delta <- as.integer(vapply(m, `[`, "", 3))

# ---- long form: one row per (rep, conditioning) ---------------------------
long <- raw %>%
  pivot_longer(c(p_path, p_union), names_to = "cond", values_to = "pval") %>%
  mutate(cond  = factor(sub("^p_", "", cond), levels = c("path", "union")),
         qlab  = factor(sprintf("q = %d", q), levels = sprintf("q = %d", sort(unique(q)))),
         dfac  = factor(delta))

# ---- plot -----------------------------------------------------------------
fig <- ggplot(long, aes(x = dfac, y = pval, fill = cond)) +
  geom_hline(yintercept = 0.05, linewidth = 0.5, linetype = "dashed",
             colour = km_ref) +
  geom_boxplot(position = position_dodge(width = 0.78), width = 0.68,
               outlier.size = 0.35, outlier.alpha = 0.35, linewidth = 0.35,
               colour = "black") +
  scale_fill_manual(values = c(path = "white", union = "#7F7F7F"),
                    labels = c(path = "Path", union = "Union"), name = NULL) +
  facet_wrap(~ qlab, nrow = 1) +
  labs(x = expression(separation~~delta), y = "Selective p-value") +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.05, 0.25, 0.5, 0.75, 1)) +
  theme_km() +
  theme(legend.key.width = unit(0.9, "cm"))

ggsave_km(fig, file.path(DIR, "..", "cw_power_box"), width = 11, height = 4.2)
cat("Wrote cw_power_box.{pdf,png} from", DIR, "\n")
