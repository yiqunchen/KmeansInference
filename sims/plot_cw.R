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
qlev <- paste0("q=", c(2, 10, 50, 100))

## Panel A: Type I QQ, FACET by q (union p-values vs Uniform) -----------------
## Single null series = the known-sigma union p-value: colour=oracle (black),
## linetype=union (dashed); q is the swept dimension -> facet, never colour.
qq <- do.call(rbind, lapply(c(2, 10, 50, 100), function(q) {
  r <- read_raw(sprintf("cw_typeI_q%d", q)); if (is.null(r) || nrow(r) < 20) return(NULL)
  pu <- sort(r$p_union[is.finite(r$p_union)])
  data.frame(q = factor(paste0("q=", q), levels = qlev),
             theo = ppoints(length(pu)), emp = pu)
}))
qq$q <- droplevels(qq$q)
pA <- ggplot(qq, aes(theo, emp)) +
  geom_abline(slope = 1, linetype = 2, colour = km_ref, alpha = 0.8) +
  geom_step(aes(colour = "oracle", linetype = "union"), linewidth = 0.9) +
  facet_wrap(~ q, nrow = 1) +
  scale_km_colour() + scale_km_linetype() +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Uniform quantile", y = "union p-value (null)",
       title = "Type I calibration") +
  theme_km()

## Panel B: power vs delta (known variance), union vs path; FACET by sigma -----
## Both curves are known-sigma -> colour=oracle (black); conditioning is the
## linetype (path solid / union dashed). sigma is the swept dimension -> facet.
pw <- s[grepl("cw_power", s$label) & s$n_valid >= 40, ]
pwk <- rbind(
  data.frame(delta = pw$delta, sig = pw$sig, rej = pw$rej_union,
             treat = "oracle", cond = "union"),
  data.frame(delta = pw$delta, sig = pw$sig, rej = pw$rej_path,
             treat = "oracle", cond = "path"))
pwk$sigf <- factor(sprintf("sigma == %s", pwk$sig))
pB <- ggplot(pwk, aes(delta, rej, colour = treat, linetype = cond)) +
  geom_hline(yintercept = 0.05, linetype = 2, colour = km_ref, alpha = 0.8) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.4, colour = km_col[["oracle"]]) +
  facet_wrap(~ sigf, nrow = 1, labeller = label_parsed) +
  scale_km_colour() + scale_km_linetype() +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = expression(paste("separation  ", delta)), y = "power",
       title = "Power vs separation (q=10)") +
  theme_km()

## Panel C: known vs unknown variance power, union vs path --------------------
## colour=treatment: known=oracle (black), unknown=studentized (blue);
## linetype=conditioning: path solid / union dashed.
c25 <- pw[pw$sig == 0.25, ]
src <- if (nrow(c25) >= 2) c25 else pw[pw$q == 10, ]   # sig=0.25 absent -> use q=10 slice
if (nrow(src) >= 2) {
  uv <- rbind(
    data.frame(delta = src$delta, rej = src$rej_union,      treat = "oracle",      cond = "union"),
    data.frame(delta = src$delta, rej = src$rej_path,       treat = "oracle",      cond = "path"),
    data.frame(delta = src$delta, rej = src$rej_union_rfib, treat = "studentized", cond = "union"),
    data.frame(delta = src$delta, rej = src$rej_path_rfib,  treat = "studentized", cond = "path"))
  uv <- uv[is.finite(uv$rej), ]
  pC <- ggplot(uv, aes(delta, rej, colour = treat, linetype = cond)) +
    geom_line(linewidth = 0.9) + geom_point(size = 2.4) +
    scale_km_colour() + scale_km_linetype() +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = expression(paste("separation  ", delta)), y = "power",
         title = "Known vs unknown variance (q=10)") +
    theme_km()
} else pC <- ggplot() + theme_void()

## ONE shared bottom legend across the three panels.
library(grid)
legC <- km_get_legend(pC + guides(colour = guide_legend(order = 1),
                                  linetype = guide_legend(order = 2)))
pA <- pA + theme(legend.position = "none")
pB <- pB + theme(legend.position = "none")
pC <- pC + theme(legend.position = "none")
body <- arrangeGrob(pA, pB, pC, ncol = 3, widths = c(1.35, 1, 1))
fig <- arrangeGrob(body, legC, ncol = 1, heights = c(10, 1.1))
fig <- arrangeGrob(fig, top = textGrob(paste0("Chen-&-Witten-faithful, n=150", tag),
                                       gp = gpar(fontface = "bold", cex = 1.0)))
ggsave_km(fig, "sims/results/cw_panels", width = 15.5, height = 4.8)
cat("Wrote sims/results/cw_panels.{pdf,png}\n")
