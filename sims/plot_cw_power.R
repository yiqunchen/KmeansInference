#!/usr/bin/env Rscript
# CW power panels (n=90, sigma=1): power vs delta stratified by q, union vs path,
# with 95% binomial CI error bars. Plus a known-vs-unknown-variance panel.
# Robust to partial data (cells with < 40 valid reps are dropped).
source("sims/house_style.R")
DIR <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "sims/results/sweep_cw"
s   <- read.csv(file.path(DIR, "summary.csv"), stringsAsFactors = FALSE)
pw  <- s[grepl("cw_power", s$label) & s$n_valid >= 40, ]
if (nrow(pw) == 0) stop("no power cells with >=40 valid reps yet")
se  <- function(p, n) sqrt(pmax(p * (1 - p), 0) / n)
done <- length(list.files(file.path(DIR, "cells"))); tag <- sprintf(" (%d/176 chunks)", done)

## --- Panel set 1: power vs delta faceted by q, union vs path (known var) ----
mk <- function(uy, py) rbind(
  data.frame(q = pw$q, delta = pw$delta, n = pw$n_valid, rej = pw[[uy]], method = "union"),
  data.frame(q = pw$q, delta = pw$delta, n = pw$n_valid, rej = pw[[py]], method = "path"))
d1 <- mk("rej_union", "rej_path")
d1$q <- factor(paste0("q = ", d1$q), levels = paste0("q = ", c(2, 10, 50)))
d1$method <- factor(d1$method, levels = c("union", "path"))
d1$lo <- pmax(0, d1$rej - 1.96 * se(d1$rej, d1$n))
d1$hi <- pmin(1, d1$rej + 1.96 * se(d1$rej, d1$n))

p1 <- ggplot(d1, aes(delta, rej, colour = method, shape = method)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.18, linewidth = 0.45, alpha = 0.8) +
  geom_line(aes(linetype = method), linewidth = 0.9) +
  geom_point(size = 2.3) +
  facet_wrap(~q) +
  scale_colour_manual(values = km_pal, labels = km_lab) +
  scale_shape_manual(values = c(union = 16, path = 17), labels = km_lab) +
  scale_linetype_manual(values = c(union = 1, path = 2), labels = km_lab) +
  scale_x_continuous(breaks = seq(1, 8, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = expression(paste("separation  ", delta)), y = "power",
       title = bquote(bold("Power by dimension q  (n=90, " * sigma * "=1, 95% CI)" * .(tag)))) +
  theme_km() +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"),
        panel.spacing = unit(1, "lines"))
ggsave_km(p1, "sims/results/cw_power", width = 11, height = 4.2)

## --- Panel set 2: known vs unknown variance (q = 10), union vs path ----------
q10 <- pw[pw$q == 10, ]
if (nrow(q10) >= 2) {
  # unknown-variance = the VALID R-fiber union/path (rfib); the phi-remap union
  # was anti-conservative and removed.
  uv <- rbind(
    data.frame(delta = q10$delta, n = q10$n_valid, rej = q10$rej_union,      method = "union", v = "known"),
    data.frame(delta = q10$delta, n = q10$n_valid, rej = q10$rej_path,       method = "path",  v = "known"),
    data.frame(delta = q10$delta, n = q10$n_valid, rej = q10$rej_union_rfib, method = "union", v = "unknown"),
    data.frame(delta = q10$delta, n = q10$n_valid, rej = q10$rej_path_rfib,  method = "path",  v = "unknown"))
  uv$method <- factor(uv$method, levels = c("union", "path"))
  uv$lo <- pmax(0, uv$rej - 1.96 * se(uv$rej, uv$n))
  uv$hi <- pmin(1, uv$rej + 1.96 * se(uv$rej, uv$n))
  p2 <- ggplot(uv, aes(delta, rej, colour = method, linetype = v, shape = v)) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.18, linewidth = 0.45, alpha = 0.8) +
    geom_line(linewidth = 0.9) + geom_point(size = 2.3) +
    scale_colour_manual(values = km_pal, labels = km_lab) +
    scale_linetype_manual(values = c(known = 1, unknown = 2), name = NULL) +
    scale_shape_manual(values = c(known = 16, unknown = 17), name = NULL) +
    scale_x_continuous(breaks = seq(1, 8, 1)) + scale_y_continuous(limits = c(0, 1)) +
    labs(x = expression(paste("separation  ", delta)), y = "power",
         title = "Known vs unknown variance (q=10)") +
    theme_km()
  ggsave_km(p2, "sims/results/cw_power_var", width = 7, height = 4.6)
}
cat("Wrote sims/results/cw_power.{pdf,png}", if (nrow(q10) >= 2) "+ cw_power_var" else "", "\n")
