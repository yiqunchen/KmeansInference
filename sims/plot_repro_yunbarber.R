#!/usr/bin/env Rscript
# FIG 1 -- "Cost of unknown sigma", APPLES-TO-APPLES (Yun-Barber's design):
# hold the truncation set FIXED (the phi-line PATH set), vary ONLY the
# sigma-handling.  Reproduces fig_power_K3_equidistant settings (q in {2,10},
# n=10/cluster, equilateral means, sigma=1, delta 0..7) in our k-means framework.
#
#   oracle (true sigma)  >=  studentized-F (same set, unknown sigma)  >=  plug-ins
#
# So unknown-sigma NEVER beats the matched oracle (the principle holds); the
# studentized-F is the best you can do without sigma and beats sigma_all; MED is
# q-dependent (craters at q=2, near-oracle at q=10).
source("sims/house_style.R")
suppressMessages({ library(ggplot2) })
load_q <- function(q) { d <- readRDS(sprintf("sims/results/repro_yunbarber_q%d.rds", q)); d$q <- q; d }
res <- rbind(load_q(2), load_q(10))
se  <- function(p, n) sqrt(pmax(p*(1-p), 0)/n)

# COLOUR = sigma-handling (km_col roles). All curves are the PATH conditioning set.
meths <- c(or_p = "oracle", uf_p = "studentized",
           md_p = "med",    sa_p = "sample")
d <- do.call(rbind, lapply(names(meths), function(col) data.frame(
  q = res$q, delta = res$delta, n = res$n_valid, rej = res[[col]], method = meths[[col]])))
d$method <- factor(d$method, levels = unname(meths))
d$qf <- factor(paste0("q = ", d$q), levels = paste0("q = ", c(2, 10)))
d$lo <- pmax(0, d$rej - 1.96*se(d$rej, d$n)); d$hi <- pmin(1, d$rej + 1.96*se(d$rej, d$n))
p <- ggplot(d, aes(delta, rej, colour = method)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, linewidth = 0.4, alpha = 0.65) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.1) +
  facet_wrap(~ qf) +
  scale_km_colour() +
  scale_x_continuous(breaks = 0:7) + scale_y_continuous(limits = c(0, 1)) +
  labs(x = expression(paste("separation  ", delta)), y = "power (path test)",
       title = "Cost of unknown variance -- same conditioning set, vary only sigma-handling") +
  theme_km() + theme(strip.background = element_blank(),
                     strip.text = element_text(face = "bold"),
                     panel.spacing = unit(0.7, "lines"))
ggsave_km(p, "sims/results/repro_yunbarber", width = 11, height = 4.8)
cat("Wrote sims/results/repro_yunbarber.{pdf,png}\n")
