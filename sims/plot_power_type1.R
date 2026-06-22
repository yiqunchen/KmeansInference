#!/usr/bin/env Rscript
# Union-vs-path Type I error + power (house style, no subtitles).
source("sims/house_style.R")

res <- readRDS("sims/results/power_type1_raw.rds")
sm  <- readRDS("sims/results/power_type1_summary.rds")
alpha <- 0.05

## Panel 1: power curves ------------------------------------------------------
## UNIFIED GRAMMAR: colour = treatment (path/union = known-sigma -> oracle black;
## naive -> orange); linetype = conditioning (path solid / union dashed; naive
## solid). Shape is not used.
pw <- data.frame(
  delta     = rep(sm$delta, 3),
  reject    = c(sm$rej_union, sm$rej_path, sm$rej_naive),
  treatment = factor(rep(c("oracle", "oracle", "naive"), each = nrow(sm)),
                     levels = c("naive", "oracle")),
  cond      = factor(rep(c("union", "path", "path"), each = nrow(sm)),
                     levels = c("path", "union")))

p1 <- ggplot(pw, aes(delta, reject, colour = treatment, linetype = cond)) +
  geom_hline(yintercept = alpha, linetype = 2, colour = km_ref, alpha = 0.6) +
  geom_line(linewidth = 1.0) +
  scale_km_colour() + scale_km_linetype() +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0.02)) +
  labs(x = expression(paste("separation  ", delta)), y = "rejection rate",
       title = "Power") +
  theme_km() + theme(legend.position = "none")

## Panel 2: null p-value ECDF vs Uniform -------------------------------------
## Only path/union here (both known-sigma -> oracle black), separated by linetype.
null0 <- res[res$delta == 0, ]
ec <- rbind(
  data.frame(p = sort(null0$p_union), F = ppoints(nrow(null0)),
             treatment = "oracle", cond = "union"),
  data.frame(p = sort(null0$p_path),  F = ppoints(nrow(null0)),
             treatment = "oracle", cond = "path"))
ec$treatment <- factor(ec$treatment, levels = c("naive", "oracle"))
ec$cond      <- factor(ec$cond, levels = c("path", "union"))

p2 <- ggplot(ec, aes(p, F, colour = treatment, linetype = cond)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = km_ref, alpha = 0.6) +
  geom_step(linewidth = 1.0) +
  scale_km_colour(drop = FALSE) + scale_km_linetype() +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "p-value", y = "empirical CDF", title = "Null calibration") +
  theme_km() + theme(legend.position = "none")

legend <- km_get_legend(p1 + theme(legend.position = "bottom"))
fig <- km_panels(list(p1, p2), legend, ncol = 2)
ggsave_km(fig, "sims/results/power_type1", width = 9, height = 4.4)
cat("Wrote sims/results/power_type1.{pdf,png}\n")
