# ---------------------------------------------------------------------------
# plot_intro_motivation.R -- intro figure following Chen-Witten (k-means) Fig 1.
# Left:   one null dataset (mu=0); k-means finds K=3 spurious clusters, centroids
#         as triangles.
# Center: QQ-plot of the NAIVE p-values over the null sweep -- grossly inflated
#         (far below the 45-degree line).
# Right:  QQ-plot of the selective (path & union) p-values -- uniform, on the line.
# ---------------------------------------------------------------------------
source("sims/house_style.R")
suppressMessages({ library(KmeansInference); library(ggplot2); library(dplyr)
                   library(gridExtra); library(grid) })

clr <- c("1" = "#CC6699", "2" = "#1F5AA6", "3" = "#E69F00")

## Left: one null dataset + spurious k-means clusters --------------------------
set.seed(5)
Xn  <- matrix(rnorm(90*2), 90, 2)
cl  <- kmeans_estimation(Xn, 3, seed = 5)$final_cluster
dfn <- data.frame(x = Xn[,1], y = Xn[,2], cl = factor(cl))
cen <- dfn %>% group_by(cl) %>% summarise(x = mean(x), y = mean(y), .groups = "drop")
pL <- ggplot(dfn, aes(x, y, colour = cl)) +
  geom_point(size = 1.9, alpha = 0.9) +
  geom_point(data = cen, aes(x, y), shape = 17, size = 4.2, colour = "black") +
  scale_colour_manual(values = clr, guide = "none") +
  labs(title = "(a)  One null dataset, K=3 clusters", x = NULL, y = NULL) +
  theme_km(12) + theme(axis.text = element_blank(), axis.ticks = element_blank())

## QQ helper + null p-values from the sweep -----------------------------------
cells <- list.files("sims/results/sweep_cw_uv/cells", pattern = "^cw_typeI_q2__c[0-9]+\\.rds$", full.names = TRUE)
r <- do.call(rbind, lapply(cells, readRDS))
qq <- function(p, lab) { p <- sort(p[is.finite(p)]); data.frame(theo = ppoints(length(p)), emp = p, lab = lab) }

pC <- ggplot(qq(r$p_naive, "naive"), aes(theo, emp)) +
  geom_abline(slope = 1, colour = km_ref, alpha = 0.7, linewidth = 0.45) +
  geom_step(colour = km_col[["naive"]], linewidth = 0.95) +
  scale_x_continuous(breaks = c(0,.5,1)) + scale_y_continuous(breaks = c(0,.5,1), limits = c(0,1)) +
  labs(title = "(b)  Naive p-values", x = "Uniform quantile", y = "sample quantile") +
  theme_km(12) + theme(aspect.ratio = 1)

seldat <- rbind(qq(r$p_path, "path"), qq(r$p_union, "union"))
seldat$lab <- factor(seldat$lab, levels = c("path", "union"))
pR <- ggplot(seldat, aes(theo, emp, colour = lab, linetype = lab)) +
  geom_abline(slope = 1, colour = km_ref, alpha = 0.7, linewidth = 0.45) +
  geom_step(linewidth = 0.95) +
  scale_colour_manual(values = c(path = km_col[["oracle"]], union = km_col[["oracle"]]),
                      labels = c(path = "Path", union = "Union"), name = NULL) +
  scale_linetype_manual(values = c(path = 1, union = 2), labels = c(path = "Path", union = "Union"), name = NULL) +
  scale_x_continuous(breaks = c(0,.5,1)) + scale_y_continuous(breaks = c(0,.5,1), limits = c(0,1)) +
  labs(title = "(c)  Selective p-values", x = "Uniform quantile", y = "sample quantile") +
  theme_km(12) + theme(aspect.ratio = 1, legend.position = c(0.72, 0.18),
                       legend.key.width = unit(0.8, "cm"), legend.background = element_blank())

g <- arrangeGrob(pL, pC, pR, ncol = 3, widths = c(1.1, 1, 1))
ggsave_km(g, "sims/results/intro_motivation", width = 12, height = 4.1)
cat("Wrote intro_motivation.{pdf,png}\n")
