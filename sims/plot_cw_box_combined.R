# ---------------------------------------------------------------------------
# plot_cw_box_combined.R -- ONE panel for BOTH calibration and power.
# x = separation delta: delta=0 is the global null (Type I), delta>0 is power.
# y = selective p-value. The fraction of each box below 0.05 is the rejection
# rate -- at delta=0 the boxes are ~Uniform (median ~0.5, calibrated); as delta
# grows they drop below 0.05 (power), and the union boxes drop faster than path.
#
# Grammar (matches cw_typeI QQ): COLOUR = variance handling (known sigma = oracle
# black, studentized-F = blue); CONDITIONING = fill, path = hollow / union =
# filled. facet_grid(variance ~ q). 0.05 + the delta=0|>0 split are grey refs.
#
# Usage: Rscript sims/plot_cw_box_combined.R [results_dir]
#   default dir = sims/results/sweep_cw_uv  (re-point to sweep_cw_uv_exact once
#   the exact power cells finish).
# ---------------------------------------------------------------------------
source("sims/house_style.R")
suppressMessages({ library(dplyr); library(tidyr) })

DIR <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "sims/results/sweep_cw_uv"
QS  <- c(2, 10, 50)   # q levels where both Type-I and power cells exist

read_cells <- function(label) {
  fs <- list.files(file.path(DIR, "cells"), pattern = paste0("^", label, "__c[0-9]+\\.rds$"),
                  full.names = TRUE)
  if (!length(fs)) return(NULL)
  do.call(rbind, lapply(fs, function(f) tryCatch(readRDS(f), error = function(e) NULL)))
}

rows <- list()
for (q in QS) {
  # delta = 0 : global-null Type-I cell ; delta >= 1 : power cells
  cells <- c(list(read_cells(sprintf("cw_typeI_q%d", q))),
             lapply(1:8, function(d) read_cells(sprintf("cw_power_q%d_d%d", q, d))))
  deltas <- c(0, 1:8)
  for (i in seq_along(cells)) {
    r <- cells[[i]]; if (is.null(r)) next
    rows[[length(rows) + 1L]] <- data.frame(
      q = q, delta = deltas[i],
      known_path = r$p_path,  known_union = r$p_union,
      stud_path  = r$p_path_rfib, stud_union = r$p_union_rfib)
  }
}
raw <- bind_rows(rows)

long <- raw %>%
  pivot_longer(c(known_path, known_union, stud_path, stud_union),
               names_to = "grp", values_to = "pval") %>%
  filter(is.finite(pval)) %>%
  mutate(
    variance = factor(ifelse(grepl("^known", grp), "Known sigma", "Studentized-F (unknown sigma)"),
                      levels = c("Known sigma", "Studentized-F (unknown sigma)")),
    cond  = factor(ifelse(grepl("union$", grp), "union", "path"), levels = c("path", "union")),
    qlab  = factor(sprintf("q = %d", q), levels = sprintf("q = %d", QS)),
    dfac  = factor(delta, levels = c(0, 1:6))) %>%
  filter(delta <= 6)

fig <- ggplot(long, aes(x = dfac, y = pval, colour = variance, fill = cond)) +
  geom_vline(xintercept = 1.5, linewidth = 0.5, linetype = "dotted", colour = km_ref) +
  geom_hline(yintercept = 0.05, linewidth = 0.45, linetype = "dashed", colour = km_ref) +
  geom_boxplot(position = position_dodge(width = 0.82), width = 0.74,
               outlier.size = 0.25, outlier.alpha = 0.3, linewidth = 0.4) +
  facet_wrap(~ qlab, nrow = 1) +
  scale_colour_manual(values = c("Known sigma" = km_col[["oracle"]],
                                 "Studentized-F (unknown sigma)" = km_col[["studentized"]]),
                      name = NULL) +
  scale_fill_manual(values = c(path = "grey90", union = "grey45"),
                    labels = c(path = "Path", union = "Union"), name = NULL) +
  labs(x = expression("separation"~~delta~~~~"("*delta*"=0: global null / Type I  |  "*delta*">0: power)"),
       y = "Selective p-value") +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.05, 0.5, 1)) +
  guides(colour = guide_legend(order = 1, override.aes = list(fill = NA, linewidth = 1.1)),
         fill   = guide_legend(order = 2, override.aes = list(colour = "black"))) +
  theme_km() +
  theme(legend.key.width = unit(1.0, "cm"))

ggsave_km(fig, file.path(DIR, "..", "cw_box_combined"), width = 12.5, height = 4.4)
cat("Wrote cw_box_combined.{pdf,png} from", DIR, "\n")
