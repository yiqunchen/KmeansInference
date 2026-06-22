#!/usr/bin/env Rscript
# Fig 1 -- Type I calibration (Chen & Witten Fig. 4 style): QQ of null p-values
# vs Uniform, facet_grid(variance ~ q). Rows = {known, unknown} variance;
# cols = q in {2,10,50,100}. UNIFIED GRAMMAR: colour = treatment / variance
# handling, linetype = conditioning (path solid / union dashed).
#   known-variance row : path & union are known-sigma selective tests -> oracle
#                        (black); naive p-value -> naive (orange).
#   unknown-variance row: path & union are STUDENTIZED (R-fiber) -> studentized
#                        (blue); naive does not appear.
# Valid tests track the 45-degree line (grey reference); naive over-rejects.
source("sims/house_style.R")
DIR <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "sims/results/sweep_cw_uv"
read_raw <- function(label) {
  fs <- list.files(file.path(DIR, "cells"), pattern = paste0("^", label, "__c"), full.names = TRUE)
  if (!length(fs)) return(NULL)
  do.call(rbind, lapply(fs, function(f) tryCatch(readRDS(f), error = function(e) NULL)))
}
qqrows <- function(p, q, cond, treatment, variance) {
  p <- sort(p[is.finite(p)]); if (length(p) < 20) return(NULL)
  data.frame(q = factor(paste0("q = ", q), levels = paste0("q = ", c(2, 10, 50, 100))),
             cond = cond, treatment = treatment, variance = variance,
             theo = ppoints(length(p)), emp = p)
}
dat <- do.call(rbind, lapply(c(2, 10, 50, 100), function(q) {
  r <- read_raw(sprintf("cw_typeI_q%d", q)); if (is.null(r)) return(NULL)
  rbind(
    # known-variance row: known-sigma selective tests (oracle) + naive
    qqrows(r$p_union,      q, "union", "oracle",      "known variance"),
    qqrows(r$p_path,       q, "path",  "oracle",      "known variance"),
    qqrows(r$p_naive,      q, "path",  "naive",       "known variance"),
    # unknown-variance row: studentized (R-fiber) tests, no naive
    qqrows(r$p_union_rfib, q, "union", "studentized", "unknown variance"),
    qqrows(r$p_path_rfib,  q, "path",  "studentized", "unknown variance"))
}))
dat$variance  <- factor(dat$variance, levels = c("known variance", "unknown variance"))
dat$cond      <- factor(dat$cond, levels = c("path", "union"))
dat$treatment <- factor(dat$treatment, levels = c("naive", "oracle", "studentized"))

p <- ggplot(dat, aes(theo, emp, colour = treatment, linetype = cond)) +
  geom_abline(slope = 1, colour = km_ref, alpha = 0.6, linewidth = 0.4) +
  geom_step(linewidth = 0.8) +
  facet_grid(variance ~ q) +
  scale_km_colour() + scale_km_linetype() +
  scale_x_continuous(breaks = c(0, 0.5, 1)) + scale_y_continuous(breaks = c(0, 0.5, 1)) +
  labs(x = "Uniform quantile", y = "selective p-value (global null)",
       title = "Type I calibration by dimension q  (n=90)") +
  theme_km() + theme(aspect.ratio = 1, panel.spacing = unit(0.6, "lines"),
                     strip.background = element_blank(),
                     strip.text = element_text(face = "bold"))
ggsave_km(p, "sims/results/cw_typeI", width = 11, height = 6)
cat("Wrote sims/results/cw_typeI.{pdf,png}\n")
