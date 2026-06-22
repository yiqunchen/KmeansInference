# Real-data example (Yun-Barber Fig 8 / Table 1 use the same data): Palmer penguins,
# 4 standardized features, k-means. For each tested cluster pair we report the
# naive p-value, the known-sigma (sigma-hat_MED plug-in) PATH and UNION p-values,
# and the studentized (unknown-sigma) PATH and UNION p-values. Headline: the
# original PATH test fails to reject obviously-real species clusters (p~0.5-0.9);
# the UNION recovers power; the studentized test is valid+powerful with no sigma.
source("sims/house_style.R")
suppressMessages({ library(KmeansInference); library(palmerpenguins); library(ggplot2); library(gridExtra) })
source("sims/kmeans_union_unknownvar.R")
K <- as.integer(Sys.getenv("K","3")); seed <- 2021
feat <- c("bill_length_mm","bill_depth_mm","flipper_length_mm","body_mass_g")
d  <- na.omit(palmerpenguins::penguins[, c("species", feat)])
X  <- scale(as.matrix(d[, feat]))
est <- kmeans_estimation(X, K, 10, seed, verbose = FALSE); cl <- est$final_cluster
s_med <- KmeansInference:::.kmeans_estimate_MED(X)
cat(sprintf("n=%d, K=%d, sigma_hat_MED=%.3f, cluster sizes: %s\n",
            nrow(X), K, s_med, paste(table(cl), collapse="/")))
prs <- combn(K, 2, simplify = FALSE)
fmt <- function(p) ifelse(p < 1e-4, sprintf("%.1e", p), sprintf("%.4f", p))
tab <- do.call(rbind, lapply(prs, function(pr) {
  c1 <- pr[1]; c2 <- pr[2]
  ft <- tryCatch(kmeans_inference_union(X, K, c1, c2, sig = s_med, seed = seed, max_paths = 3000), error=function(e) NULL)
  rf <- tryCatch(kmeans_union_unknownvar(X, K, c1, c2, seed = seed), error=function(e) NULL)
  if (is.null(ft) || is.null(rf)) return(NULL)
  data.frame(pair = sprintf("%d vs %d", c1, c2), n1 = sum(cl==c1), n2 = sum(cl==c2),
             naive = ft$p_naive, path_known = ft$pval_path, union_known = ft$pval_union,
             path_studF = rf$p_path, union_studF = rf$p_union)
}))
out <- tab; for (j in 4:8) out[[j]] <- fmt(tab[[j]])
cat("\n=== Penguins cluster-pair selective p-values ===\n"); print(out, row.names = FALSE)
write.csv(tab, sprintf("sims/results/penguins_pvalues_k%d.csv", K), row.names = FALSE)
cat(sprintf("Wrote sims/results/penguins_pvalues_k%d.csv  (figure: sims/plot_penguins.R)\n", K))
