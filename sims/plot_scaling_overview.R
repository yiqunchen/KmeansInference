#!/usr/bin/env Rscript
# Cluster-sizing overview from the existing CW sweep (n=90, k=3): per-rep TOTAL
# runtime and Lloyd-path count vs dimension q, one line per separation delta.
# (Total runtime = all methods in one rep: known union/path + R-fiber + plug-ins.)
# Takeaway: cost tracks the number of Lloyd paths explored, which is LARGEST at
# small q and small delta -- not at large q.  Use mean_n_paths as the cost driver.
#
# GRAMMAR EXCEPTION (documented): this is a pure COST SWEEP with NO inference
# treatment and NO path/union conditioning series -- only the swept separation
# delta. Per the unified grammar, a sequential (viridis) colour for the swept
# dimension delta is the allowed exception here; x = the other swept dimension q.
source("sims/house_style.R")
suppressMessages({ library(ggplot2); library(gridExtra); library(grid) })
DIR <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "sims/results/sweep_cw_uv"
s  <- read.csv(file.path(DIR, "summary.csv"), stringsAsFactors = FALSE)
pw <- s[grepl("cw_power", s$label) & s$n_valid >= 40, ]
pw$dl <- factor(pw$delta)

base <- function(y, ylab, ttl) ggplot(pw, aes(factor(q), .data[[y]], colour = dl, group = dl)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  scale_colour_viridis_d(name = expression(delta), option = "C", end = 0.92) +
  labs(x = "dimension  q", y = ylab, title = ttl) + theme_km()
p1 <- base("median_runtime_s", "median runtime / rep (s)", "Compute time (all methods, 1 rep)")
p2 <- base("mean_n_paths",     "mean # Lloyd paths",       "Path-exploration complexity")
g <- arrangeGrob(p1, p2, ncol = 2,
                 top = textGrob("CW sweep cost overview (n=90, k=3) -- cost tracks path count, peaks at small q/small delta",
                                gp = gpar(fontface = "bold", fontsize = 12)))
ggsave_km(g, file.path(dirname(DIR), "scaling_overview"), width = 11, height = 4.6)
cat("Wrote scaling_overview.{pdf,png} (in", dirname(DIR), ")\n")
