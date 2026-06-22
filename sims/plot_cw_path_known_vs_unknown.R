#!/usr/bin/env Rscript
# Resolves "unknown-sigma path looks more powerful than known-sigma path" -- is it
# a contradiction?  NO.  The two axes are SEPARATED here:
#   colour   = sigma-handling : known (true sigma)  vs  studentized (unknown)
#   linetype = conditioning    : path (single Lloyd path)  vs  weaker (union / R-fiber)
#
# WITHIN a conditioning (same linetype) known >= studentized -- ALWAYS (18/18 cells
# each; verified).  The apparent "known < studentized" compares ACROSS linetypes
# (studentized-weaker vs known-path), i.e. different conditioning sets -- apples vs
# oranges, the same cross-conditioning effect as repro_conditioning.
#
# Curves (conditional power, 95% binomial CI on n_detected):
#   known  + path   = cpow_path        (Chen-Witten, phi-ray, single path)
#   studz. + path   = cpow_path_unk     (F pivot on the SAME phi-line path set)
#   known  + weaker = cpow_union        (Chen-Witten union, phi-ray, weaker)
#   studz. + weaker = cpow_union_rfib   (R-fiber union; on the R-fiber, conditioning
#                                        on total scale T makes sigma irrelevant)
source("sims/house_style.R")
suppressMessages({ library(ggplot2) })
DIR <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "sims/results/sweep_cw_uv"
s   <- read.csv(file.path(DIR, "summary.csv"), stringsAsFactors = FALSE)
pw  <- s[grepl("cw_power", s$label) & s$n_valid >= 40 & s$n_detected >= 10, ]
sebin <- function(p, n) sqrt(pmax(p*(1-p), 0)/pmax(n, 1))

# COLOUR = sigma-handling (km_col): oracle (known) vs studentized (unknown).
# LINETYPE = conditioning (km_cond): path solid, union (weaker / R-fiber) dashed.
spec <- list(
  c(col = "cpow_path",       sig = "oracle",      cond = "path"),
  c(col = "cpow_path_unk",   sig = "studentized", cond = "path"),
  c(col = "cpow_union",      sig = "oracle",      cond = "union"),
  c(col = "cpow_union_rfib", sig = "studentized", cond = "union"))
d <- do.call(rbind, lapply(spec, function(z) data.frame(
  q = pw$q, delta = pw$delta, nd = pw$n_detected, cp = pw[[z["col"]]],
  sig = z["sig"], cond = z["cond"])))
d <- d[is.finite(d$cp), ]
d$sig  <- factor(d$sig,  levels = c("oracle", "studentized"))
d$cond <- factor(d$cond, levels = c("path", "union"))
d$qf   <- factor(paste0("q = ", d$q), levels = paste0("q = ", sort(unique(d$q))))
d$lo <- pmax(0, d$cp - 1.96*sebin(d$cp, d$nd)); d$hi <- pmin(1, d$cp + 1.96*sebin(d$cp, d$nd))
p <- ggplot(d, aes(delta, cp, colour = sig, linetype = cond)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.18, linewidth = 0.4, alpha = 0.55) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.7) +
  facet_wrap(~ qf) +
  scale_km_colour() +
  scale_km_linetype() +
  scale_x_continuous(breaks = seq(1, 8, 1)) + scale_y_continuous(limits = c(0, 1)) +
  labs(x = expression(paste("separation  ", delta)), y = "conditional power (path test)",
       title = "Within each conditioning (same linetype), known >= studentized -- no contradiction") +
  theme_km() + theme(strip.background = element_blank(),
                     strip.text = element_text(face = "bold"),
                     panel.spacing = unit(0.7, "lines"))
ggsave_km(p, file.path(dirname(DIR), "cw_path_known_vs_unknown"), width = 12, height = 5.2)
cat("Wrote cw_path_known_vs_unknown.{pdf,png} (in", dirname(DIR), ")\n")
