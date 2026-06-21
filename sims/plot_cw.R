#!/usr/bin/env Rscript
# Chen-&-Witten-faithful (n=150) panels for the CW sweep -- house style.
# Robust to PARTIAL data (sweep may still be running): cells/curves with too
# few valid reps are dropped, and the title is marked accordingly.
source("sims/house_style.R")
DIR <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "sims/results/sweep_cw"
s   <- read.csv(file.path(DIR, "summary.csv"), stringsAsFactors = FALSE)
done <- length(list.files(file.path(DIR, "cells")))
tag  <- sprintf(" (partial: %d/140 chunks)", done)

read_raw <- function(label) {
  fs <- list.files(file.path(DIR, "cells"),
                   pattern = paste0("^", label, "__c"), full.names = TRUE)
  if (!length(fs)) return(NULL)
  do.call(rbind, lapply(fs, function(f) tryCatch(readRDS(f), error = function(e) NULL)))
}
qcol <- c("q=2" = "#1F5AA6", "q=10" = "#C98A2E", "q=50" = "#2B7A78", "q=100" = "#9A4D8E")

## Panel A: Type I QQ across q (union p-values vs Uniform) --------------------
qq <- do.call(rbind, lapply(c(2, 10, 50, 100), function(q) {
  r <- read_raw(sprintf("cw_typeI_q%d", q)); if (is.null(r) || nrow(r) < 20) return(NULL)
  pu <- sort(r$p_union[is.finite(r$p_union)])
  data.frame(q = factor(paste0("q=", q), levels = names(qcol)),
             theo = ppoints(length(pu)), emp = pu)
}))
pA <- ggplot(qq, aes(theo, emp, colour = q)) +
  geom_abline(slope = 1, linetype = 2, colour = km_ref, alpha = 0.6) +
  geom_step(linewidth = 0.9) +
  scale_colour_manual(values = qcol) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Uniform quantile", y = "union p-value (null)",
       title = "Type I calibration") +
  theme_km() + theme(legend.position = c(0.02, 0.98),
                     legend.justification = c(0, 1))

## Panel B: power vs delta by sigma (known variance), union vs path -----------
pw <- s[grepl("cw_power", s$label) & s$n_valid >= 40, ]
scol <- c("0.25" = "#1F5AA6", "0.5" = "#C98A2E", "1" = "#B55D4C")
mk_long <- function(uy, py, lab) rbind(
  data.frame(delta = pw$delta, sig = factor(pw$sig), rej = pw[[uy]], method = "union", panel = lab),
  data.frame(delta = pw$delta, sig = factor(pw$sig), rej = pw[[py]], method = "path",  panel = lab))
pwk <- mk_long("rej_union", "rej_path", "known")
pB <- ggplot(pwk, aes(delta, rej, colour = sig, linetype = method, shape = method)) +
  geom_hline(yintercept = 0.05, linetype = 3, colour = km_ref, alpha = 0.6) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.4) +
  scale_colour_manual(values = scol, name = expression(sigma)) +
  scale_linetype_manual(values = c(union = 1, path = 2)) +
  scale_shape_manual(values = c(union = 16, path = 17)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = expression(paste("separation  ", delta)), y = "power",
       title = "Power vs separation (q=10)") +
  theme_km()

## Panel C: known vs unknown variance power (sigma = 0.25), union vs path -----
c25 <- pw[pw$sig == 0.25, ]
if (nrow(c25) >= 2) {
  uv <- rbind(
    data.frame(delta = c25$delta, rej = c25$rej_union,     method = "union", v = "known"),
    data.frame(delta = c25$delta, rej = c25$rej_path,      method = "path",  v = "known"),
    data.frame(delta = c25$delta, rej = c25$rej_union_unk, method = "union", v = "unknown"),
    data.frame(delta = c25$delta, rej = c25$rej_path_unk,  method = "path",  v = "unknown"))
  uv$method <- factor(uv$method, levels = c("union", "path"))
  pC <- ggplot(uv, aes(delta, rej, colour = method, linetype = v, shape = v)) +
    geom_line(linewidth = 0.9) + geom_point(size = 2.4) +
    scale_colour_manual(values = km_pal, labels = km_lab) +
    scale_linetype_manual(values = c(known = 1, unknown = 2), name = NULL) +
    scale_shape_manual(values = c(known = 16, unknown = 17), name = NULL) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = expression(paste("separation  ", delta)), y = "power",
         title = expression(paste("Known vs unknown variance (", sigma, "=0.25)"))) +
    theme_km()
} else pC <- ggplot() + theme_void()

library(grid)
fig <- arrangeGrob(pA, pB, pC, ncol = 3)
fig <- arrangeGrob(fig, top = textGrob(paste0("Chen-&-Witten-faithful, n=150", tag),
                                       gp = gpar(fontface = "bold", cex = 1.0)))
ggsave_km(fig, "sims/results/cw_panels", width = 14, height = 4.6)
cat("Wrote sims/results/cw_panels.{pdf,png}\n")
