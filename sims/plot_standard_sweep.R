#!/usr/bin/env Rscript
# Standard-sweep summary: where the union gains, and the extension settings
# (k-sweep, general covariance, unbalanced) -- house style.
source("sims/house_style.R")
s <- read.csv("sims/results/sweep_standard/summary.csv", stringsAsFactors = FALSE)

## Panel 1: power gain (union - path) vs delta, by q (sigma=1) ----------------
g <- s[grepl("^power_q.*_s1_d", s$label), ]
g$q <- factor(paste0("q=", g$q))
p1 <- ggplot(g, aes(delta, power_gain, colour = q, shape = q)) +
  geom_hline(yintercept = 0, linetype = 2, colour = km_ref, alpha = 0.6) +
  geom_line(linewidth = 1.0) + geom_point(size = 2.8) +
  scale_colour_manual(values = c("q=2" = "#1F5AA6", "q=10" = "#C98A2E")) +
  scale_shape_manual(values = c("q=2" = 16, "q=10" = 17)) +
  labs(x = expression(paste("separation  ", delta)),
       y = "power gain  (union - path)", title = "Where the union gains") +
  theme_km()

## Panel 2: union vs path power across extension settings (delta=5) -----------
ext <- s[s$label %in% c("k2", "k4", "k5", "unbal_d5",
                        "gencov_ar0.5_d5", "gencov_ar0.9_d5"), ]
nice <- c(k2 = "k=2", k4 = "k=4", k5 = "k=5", unbal_d5 = "unbalanced",
          `gencov_ar0.5_d5` = "AR(0.5)", `gencov_ar0.9_d5` = "AR(0.9)")
ext$name <- factor(nice[ext$label], levels = rev(nice))
pe <- rbind(
  data.frame(name = ext$name, power = ext$rej_path,  method = "path"),
  data.frame(name = ext$name, power = ext$rej_union, method = "union"))
pe$method <- factor(pe$method, levels = c("union", "path"))
p2 <- ggplot(pe, aes(power, name, colour = method)) +
  geom_line(aes(group = name), colour = "grey70", linewidth = 1.0) +
  geom_point(size = 3.4) +
  scale_colour_manual(values = km_pal, labels = km_lab) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(x = "power", y = NULL, title = expression(paste("Extensions (", delta, " = 5)"))) +
  theme_km() + theme(legend.position = "none")

legend <- km_get_legend(p2 + theme(legend.position = "bottom"))
# panel 1 has its own (q) legend; keep it, and add the method legend under panel 2
library(grid)
fig <- arrangeGrob(p1, arrangeGrob(p2, legend, ncol = 1, heights = c(10, 1.4)),
                   ncol = 2)
ggsave_km(fig, "sims/results/standard_sweep", width = 10, height = 4.6)
cat("Wrote sims/results/standard_sweep.{pdf,png}\n")
