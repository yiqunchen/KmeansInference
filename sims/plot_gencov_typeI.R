#!/usr/bin/env Rscript
# Type I calibration under NON-ISOTROPIC (AR(1)) covariance, supplying the correct
# SigInv (Yun-Barber Fig 7 spirit; ours is a stronger AR test than their diagonal).
# Global null (delta=0), k=3, q=2, AR(0.5) and AR(0.9). Known-sigma union/path stay
# on the 45-degree line; the naive (non-selective) p-value over-rejects.
# Grammar: colour = treatment (naive orange / known-oracle black); linetype =
# conditioning (path solid, union dashed); grey = 45-degree reference.
source("sims/house_style.R"); suppressMessages(library(ggplot2))
read_raw <- function(lab) {
  fs <- list.files("sims/results/sweep_standard/cells", paste0("^", lab, "__c"), full.names = TRUE)
  if (!length(fs)) return(NULL)
  do.call(rbind, lapply(fs, function(f) tryCatch(readRDS(f), error = function(e) NULL)))
}
qqrows <- function(p, cov, treat, cond) {
  p <- sort(p[is.finite(p)]); if (length(p) < 50) return(NULL)
  data.frame(cov = cov, treat = treat, cond = cond, theo = ppoints(length(p)), emp = p)
}
covs <- c(gencov_ar0.5_d0 = "AR(0.5) covariance", gencov_ar0.9_d0 = "AR(0.9) covariance")
dat <- do.call(rbind, lapply(names(covs), function(lab) {
  r <- read_raw(lab); if (is.null(r)) return(NULL)
  rbind(qqrows(r$p_path,  covs[[lab]], "oracle", "path"),
        qqrows(r$p_union, covs[[lab]], "oracle", "union"),
        qqrows(r$p_naive, covs[[lab]], "naive",  "path"))
}))
dat$cov   <- factor(dat$cov, levels = unname(covs))
dat$treat <- factor(dat$treat, levels = c("naive", "oracle"))
dat$cond  <- factor(dat$cond, levels = c("path", "union"))
rate <- aggregate(emp ~ cov + treat + cond, dat, function(p) mean(p <= 0.05))
cat("Empirical Type I (p<=0.05):\n"); print(rate, row.names = FALSE, digits = 3)
p <- ggplot(dat, aes(theo, emp, colour = treat, linetype = cond)) +
  geom_abline(slope = 1, colour = km_ref, alpha = 0.55, linewidth = 0.4) +
  geom_line(linewidth = 0.85) +
  facet_wrap(~ cov) +
  scale_km_colour() + scale_km_linetype() +
  scale_x_continuous(breaks = c(0, 0.5, 1)) + scale_y_continuous(breaks = c(0, 0.5, 1)) +
  labs(x = "Uniform quantile", y = "selective p-value (global null)",
       title = "Type I under non-isotropic AR covariance (correct SigInv)") +
  theme_km() + theme(aspect.ratio = 1, panel.spacing = unit(0.6, "lines"))
ggsave_km(p, "sims/results/gencov_typeI", width = 9, height = 4.6)
cat("Wrote sims/results/gencov_typeI.{pdf,png}\n")
