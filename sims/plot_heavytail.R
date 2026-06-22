#!/usr/bin/env Rscript
# Heavy-tailed robustness (Yun-Barber Figs 5-6): global-null Type I QQ when noise is
# Gaussian / t5 / t10 (unit variance). facet_grid(noise ~ q). Grammar: colour =
# treatment (naive orange / known-oracle black / studentized blue), linetype =
# conditioning (path solid, union dashed), grey = 45-degree reference.
source("sims/house_style.R"); suppressMessages(library(ggplot2))
r <- readRDS("sims/results/heavytail_typeI.rds")
qqrows <- function(p, noise, q, treat, cond) {
  p <- sort(p[is.finite(p)]); if (length(p) < 50) return(NULL)
  data.frame(noise = noise, q = q, treat = treat, cond = cond, theo = ppoints(length(p)), emp = p)
}
spec <- list(c("p_path","oracle","path"), c("p_union","oracle","union"),
             c("p_path_rfib","studentized","path"), c("p_union_rfib","studentized","union"),
             c("p_naive","naive","path"))
dat <- do.call(rbind, lapply(split(r, list(r$noise, r$q), drop = TRUE), function(g) {
  do.call(rbind, lapply(spec, function(s) qqrows(g[[s[1]]], g$noise[1], g$q[1], s[2], s[3]))) }))
dat$noise <- factor(dat$noise, levels = c("gaussian","t5","t10"),
                    labels = c("Gaussian","t5 noise","t10 noise"))
dat$qf    <- factor(paste0("q = ", dat$q), levels = paste0("q = ", c(2,10)))
dat$treat <- factor(dat$treat, levels = c("naive","oracle","studentized"))
dat$cond  <- factor(dat$cond, levels = c("path","union"))
p <- ggplot(dat, aes(theo, emp, colour = treat, linetype = cond)) +
  geom_abline(slope = 1, colour = km_ref, alpha = 0.55, linewidth = 0.4) +
  geom_line(linewidth = 0.8) +
  facet_grid(noise ~ qf) +
  scale_km_colour() + scale_km_linetype() +
  scale_x_continuous(breaks = c(0, 0.5, 1)) + scale_y_continuous(breaks = c(0, 0.5, 1)) +
  labs(x = "Uniform quantile", y = "selective p-value (global null)",
       title = "Type I under heavy-tailed noise (t5, t10): selective tests calibrated, naive invalid") +
  theme_km() + theme(aspect.ratio = 1, panel.spacing = unit(0.5, "lines"))
ggsave_km(p, "sims/results/heavytail_typeI", width = 9, height = 9)
cat("Wrote sims/results/heavytail_typeI.{pdf,png}\n")
