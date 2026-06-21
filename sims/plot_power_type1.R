#!/usr/bin/env Rscript
# Union-vs-path Type I error + power (house style, no subtitles).
source("sims/house_style.R")

res <- readRDS("sims/results/power_type1_raw.rds")
sm  <- readRDS("sims/results/power_type1_summary.rds")
alpha <- 0.05

## Panel 1: power curves ------------------------------------------------------
pw <- data.frame(
  delta  = rep(sm$delta, 3),
  reject = c(sm$rej_union, sm$rej_path, sm$rej_naive),
  method = factor(rep(c("union", "path", "naive"), each = nrow(sm)),
                  levels = c("union", "path", "naive")))

p1 <- ggplot(pw, aes(delta, reject, colour = method, shape = method)) +
  geom_hline(yintercept = alpha, linetype = 2, colour = km_ref, alpha = 0.6) +
  geom_line(linewidth = 1.0) + geom_point(size = 2.8) +
  scale_colour_manual(values = km_pal, labels = km_lab) +
  scale_shape_manual(values = c(union = 16, path = 16, naive = 4), labels = km_lab) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0.02)) +
  labs(x = expression(paste("separation  ", delta)), y = "rejection rate",
       title = "Power") +
  theme_km() + theme(legend.position = "none")

## Panel 2: null p-value ECDF vs Uniform -------------------------------------
null0 <- res[res$delta == 0, ]
ec <- rbind(
  data.frame(p = sort(null0$p_union), F = ppoints(nrow(null0)), method = "union"),
  data.frame(p = sort(null0$p_path),  F = ppoints(nrow(null0)), method = "path"))
ec$method <- factor(ec$method, levels = c("union", "path", "naive"))

p2 <- ggplot(ec, aes(p, F, colour = method)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = km_ref, alpha = 0.6) +
  geom_step(linewidth = 1.0) +
  scale_colour_manual(values = km_pal, labels = km_lab, drop = FALSE) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "p-value", y = "empirical CDF", title = "Null calibration") +
  theme_km() + theme(legend.position = "none")

legend <- km_get_legend(p1 + theme(legend.position = "bottom"))
fig <- km_panels(list(p1, p2), legend, ncol = 2)
ggsave_km(fig, "sims/results/power_type1", width = 9, height = 4.4)
cat("Wrote sims/results/power_type1.{pdf,png}\n")
