#!/usr/bin/env Rscript
# FIG 2 -- Union vs path is a CONDITIONING effect, orthogonal to the variance axis.
# Same settings (q in {2,10}, n=10/cluster, equilateral, sigma=1, delta 0..7).
#
# colour = sigma knowledge: known (true sigma, phi-ray) vs unknown (R-fiber F).
# linetype/shape = conditioning set: path (single clustering) vs union (Lloyd paths).
#
# Read-offs:
#  * KNOWN sigma: union (grey solid) >> path (grey dashed)  -> the union power gain.
#  * UNKNOWN sigma (R-fiber): union (blue solid) ~= path (blue dashed)  -> the gain
#    is "studentized away".
#  * LADDER: the R-fiber unknown-sigma curves sit BETWEEN known-path and known-union
#    -- weaker conditioning lifts the R-fiber test ABOVE the known-sigma PATH oracle
#    NOT because unknown beats known, but because it conditions on less. (This is why
#    comparing R-fiber-F to the phi-ray oracle earlier looked paradoxical.)
source("sims/house_style.R")
suppressMessages({ library(ggplot2) })
load_q <- function(q) { d <- readRDS(sprintf("sims/results/repro_yunbarber_q%d.rds", q)); d$q <- q; d }
res <- rbind(load_q(2), load_q(10))
se  <- function(p, n) sqrt(pmax(p*(1-p), 0)/n)

# COLOUR = sigma-handling (km_col): oracle (known) vs studentized (unknown, R-fiber).
# LINETYPE = conditioning (km_cond): path solid, union dashed.
spec <- list(
  c(col = "or_p", sig = "oracle",      cond = "path"),
  c(col = "or_u", sig = "oracle",      cond = "union"),
  c(col = "rf_p", sig = "studentized", cond = "path"),
  c(col = "rf_u", sig = "studentized", cond = "union"))
d <- do.call(rbind, lapply(spec, function(s) data.frame(
  q = res$q, delta = res$delta, n = res$n_valid, rej = res[[s["col"]]],
  sig = s["sig"], cond = s["cond"])))
d$sig  <- factor(d$sig,  levels = c("oracle", "studentized"))
d$cond <- factor(d$cond, levels = c("path", "union"))
d$qf   <- factor(paste0("q = ", d$q), levels = paste0("q = ", c(2, 10)))
d$lo <- pmax(0, d$rej - 1.96*se(d$rej, d$n)); d$hi <- pmin(1, d$rej + 1.96*se(d$rej, d$n))
p <- ggplot(d, aes(delta, rej, colour = sig, linetype = cond)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.1) +
  facet_wrap(~ qf) +
  scale_km_colour() +
  scale_km_linetype() +
  scale_x_continuous(breaks = 0:7) + scale_y_continuous(limits = c(0, 1)) +
  labs(x = expression(paste("separation  ", delta)), y = "power",
       title = "Union vs path: a conditioning gain (known sigma) that is studentized away (unknown sigma)") +
  theme_km() + theme(strip.background = element_blank(),
                     strip.text = element_text(face = "bold"),
                     panel.spacing = unit(0.7, "lines"))
ggsave_km(p, "sims/results/repro_conditioning", width = 11, height = 4.8)
cat("Wrote sims/results/repro_conditioning.{pdf,png}\n")
