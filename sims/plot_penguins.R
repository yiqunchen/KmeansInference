#!/usr/bin/env Rscript
# Penguins real-data figure (reads the p-value CSV from real_data_penguins.R; only
# re-fits the fast k-means for the scatter). Panel A: PCA scatter coloured by
# k-means cluster, shaped by true species. Panel B: -log10 selective p-value by
# cluster pair -- colour = treatment (sigma_MED plug-in green vs studentized blue),
# shape = conditioning (path open / union filled); naive omitted (off-scale invalid).
source("sims/house_style.R")
suppressMessages({ library(KmeansInference); library(palmerpenguins); library(ggplot2); library(gridExtra) })
K <- as.integer(Sys.getenv("K","3")); seed <- 2021
feat <- c("bill_length_mm","bill_depth_mm","flipper_length_mm","body_mass_g")
d <- na.omit(palmerpenguins::penguins[, c("species", feat)])
X <- scale(as.matrix(d[, feat]))
cl <- kmeans_estimation(X, K, 10, seed, verbose = FALSE)$final_cluster
tab <- read.csv(sprintf("sims/results/penguins_pvalues_k%d.csv", K), stringsAsFactors = FALSE)

## Panel A: PCA scatter (cluster = qualitative data palette = documented exception)
pcs <- prcomp(X)$x[, 1:2]
clpal <- c("#882255","#44AA99","#DDCC77","#332288","#AA4499")[seq_len(K)]; names(clpal) <- as.character(seq_len(K))
fdf <- data.frame(PC1 = pcs[,1], PC2 = pcs[,2], cluster = factor(cl), species = d$species)
pA <- ggplot(fdf, aes(PC1, PC2, colour = cluster, shape = species)) +
  geom_point(size = 1.9, alpha = 0.85) +
  scale_colour_manual(values = clpal, name = "k-means cluster") +
  scale_shape_discrete(name = "true species") +
  labs(title = sprintf("Palmer penguins: k-means (K=%d)", K)) +
  theme_km() + theme(legend.position = "right",
                     legend.title = element_text(face = "bold", size = 9), legend.text = element_text(size = 8))

## Panel B: -log10 p by pair (x = pair discrete, y = nlp continuous -> standard dodge)
nlp <- function(p) -log10(pmax(p, 1e-300))
mlt <- do.call(rbind, lapply(seq_len(nrow(tab)), function(i) data.frame(
  pair  = tab$pair[i],
  treat = c("med","med","studentized","studentized"),
  cond  = c("path","union","path","union"),
  nlp   = nlp(c(tab$path_known[i], tab$union_known[i], tab$path_studF[i], tab$union_studF[i])))))
mlt$treat <- factor(mlt$treat, levels = c("med","studentized"))
mlt$cond  <- factor(mlt$cond,  levels = c("path","union"))
pB <- ggplot(mlt, aes(pair, nlp, colour = treat, shape = cond, group = interaction(treat, cond))) +
  geom_hline(yintercept = nlp(0.05), linetype = 2, colour = km_ref, alpha = 0.8) +
  geom_point(size = 3, position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = c(med = km_col[["med"]], studentized = km_col[["studentized"]]),
                      labels = c(med = "sigma-MED plug-in", studentized = "Studentized-F"), name = NULL) +
  scale_shape_manual(values = c(path = 1, union = 16), labels = km_cond_lab, name = NULL) +
  annotate("text", x = 0.6, y = nlp(0.05) + 0.5, label = "p = 0.05", size = 2.7, colour = km_ref, hjust = 0) +
  labs(x = "cluster pair", y = expression(paste(-log[10], "  selective p-value")),
       title = "Studentized rejects all 3; known-sigma union only 2v3") +
  theme_km()
g <- arrangeGrob(pA, pB, ncol = 2, widths = c(1.25, 1))
ggsave_km(g, "sims/results/penguins", width = 13.5, height = 4.6)
cat("Wrote sims/results/penguins.{pdf,png}\n")
