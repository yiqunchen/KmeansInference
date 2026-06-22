#!/usr/bin/env Rscript
# Union vs path across STRUCTURAL settings (number of clusters K, cluster balance,
# covariance) at delta=5, sigma=1, n per cluster=20 -- a Cleveland dot-plot. The
# path->union gap IS the known-sigma power gain; it WIDENS with K (k=2..5). (The
# power-vs-delta view lives in cw_power; this figure carries only the extensions,
# and is the single home for the K-sweep.)
# GRAMMAR EXCEPTION (documented in house_style.R): the dot-plot uses point SHAPE
# (path = open 1, union = filled 16) because linetype cannot separate two points on
# one row; single oracle-black colour (all known-sigma).
source("sims/house_style.R"); suppressMessages(library(grid))
s <- read.csv("sims/results/sweep_standard/summary.csv", stringsAsFactors = FALSE)
# K=3 is the matching power_q2_s1_d5 cell (q=2, delta=5, sigma=1); relabel to k3
k3  <- transform(s[s$label == "power_q2_s1_d5", ], label = "k3")
ext <- rbind(s[s$label %in% c("k2","k4","k5","unbal_d5","gencov_ar0.5_d5","gencov_ar0.9_d5"), ], k3)
nice <- c(k2="k=2 clusters", k3="k=3 clusters", k4="k=4 clusters", k5="k=5 clusters",
          unbal_d5="unbalanced sizes", `gencov_ar0.5_d5`="AR(0.5) covariance",
          `gencov_ar0.9_d5`="AR(0.9) covariance")
ext$name <- factor(nice[ext$label], levels = rev(nice[c("k2","k3","k4","k5","unbal_d5","gencov_ar0.5_d5","gencov_ar0.9_d5")]))
pe <- rbind(data.frame(name = ext$name, power = ext$rej_path,  cond = "path"),
            data.frame(name = ext$name, power = ext$rej_union, cond = "union"))
p <- ggplot(pe, aes(power, name)) +
  geom_line(aes(group = name), colour = "grey70", linewidth = 1.0) +
  geom_point(aes(shape = cond), colour = km_col[["oracle"]], size = 3.4) +
  scale_shape_manual(values = c(path = 1, union = 16), labels = km_cond_lab, name = NULL) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(x = "power (known sigma)", y = NULL,
       title = expression(paste("Union vs path across settings  (", delta, "=5, n/cluster=20): gap = the union gain"))) +
  theme_km() + theme(legend.position = "bottom")
ggsave_km(p, "sims/results/standard_sweep", width = 7, height = 4.4)
cat("Wrote sims/results/standard_sweep.{pdf,png}\n")
