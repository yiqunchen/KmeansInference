#!/usr/bin/env Rscript
# Type I calibration QQ, stratified by q, for union AND path, under known AND
# unknown variance. Encoding (one channel = one variable): colour = q,
# linetype = method (union solid / path dashed); facet = variance.
source("sims/house_style.R")
DIR <- "sims/results/sweep_cw"
read_raw <- function(label) {
  fs <- list.files(file.path(DIR, "cells"),
                   pattern = paste0("^", label, "__c"), full.names = TRUE)
  if (!length(fs)) return(NULL)
  do.call(rbind, lapply(fs, function(f) tryCatch(readRDS(f), error = function(e) NULL)))
}
qcol <- c("q=2" = "#1F5AA6", "q=10" = "#C98A2E", "q=50" = "#2B7A78", "q=100" = "#9A4D8E")

qqrows <- function(p, q, method, variance) {
  p <- sort(p[is.finite(p)]); if (length(p) < 20) return(NULL)
  data.frame(q = factor(paste0("q=", q), levels = names(qcol)),
             method = method, variance = variance,
             theo = ppoints(length(p)), emp = p)
}
dat <- do.call(rbind, lapply(c(2, 10, 50, 100), function(q) {
  r <- read_raw(sprintf("cw_typeI_q%d", q)); if (is.null(r)) return(NULL)
  rbind(qqrows(r$p_union,         q, "union", "known variance"),
        qqrows(r$p_path,          q, "path",  "known variance"),
        qqrows(r$p_union_unknown, q, "union", "unknown variance"),
        qqrows(r$p_path_unknown,  q, "path",  "unknown variance"))
}))
dat$variance <- factor(dat$variance, levels = c("known variance", "unknown variance"))
dat$method   <- factor(dat$method, levels = c("union", "path"))

n_reps <- max(table(interaction(dat$q, dat$method, dat$variance))) # for the title

p <- ggplot(dat, aes(theo, emp, colour = q, linetype = method)) +
  geom_abline(slope = 1, colour = km_ref, alpha = 0.6, linewidth = 0.4) +
  geom_step(linewidth = 0.8) +
  facet_wrap(~variance) +
  scale_colour_manual(values = qcol, name = NULL) +
  scale_linetype_manual(values = c(union = 1, path = 2), name = NULL) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Uniform quantile", y = "selective p-value (global null)",
       title = "Type I calibration by dimension q (n=150)") +
  theme_km() + theme(legend.position = "right",
                     panel.spacing = unit(1, "lines"),
                     strip.background = element_blank(),
                     strip.text = element_text(face = "bold"))

ggsave_km(p, "sims/results/cw_typeI", width = 10, height = 5)
cat("Wrote sims/results/cw_typeI.{pdf,png}\n")
