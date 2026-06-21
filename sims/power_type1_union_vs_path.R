#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Type I error + power comparison: union test vs. original (path) test.
#
# Logic of the check
# ------------------
# kmeans_inference_union() returns, from a SINGLE fit, three p-values for the
# same contrast on the same data:
#   * pval_union : the new union-over-Lloyd-paths test
#   * pval_path  : the original kmeans_inference() p-value (verified identical
#                  by tests/testthat/test-kmeans-union.R)
#   * p_naive    : the invalid z-test that ignores selection (reference)
#
# Predictions:
#   delta = 0 (global null, all cluster means equal) -> H0 true for every pair
#       => p_union and p_path ~ Uniform(0,1); reject rate ~ alpha.   (Type I)
#   delta > 0 (separated clusters)
#       => p_union should reject MORE often than p_path (more power),
#          both << p_naive's over-rejection.                          (Power)
#
# Run: Rscript sims/power_type1_union_vs_path.R [nreps] [max_paths]
# ---------------------------------------------------------------------------

suppressMessages(library(KmeansInference))
suppressMessages(library(parallel))

args      <- commandArgs(trailingOnly = TRUE)
nreps     <- if (length(args) >= 1) as.integer(args[1]) else 200L
max_paths <- if (length(args) >= 2) as.integer(args[2]) else 100L
ncores    <- if (length(args) >= 3) as.integer(args[3]) else max(1L, detectCores() - 1L)

n_per   <- 10L          # observations per true cluster
k       <- 3L
q       <- 2L
sig     <- 1
iter.max<- 10L
deltas  <- c(0, 1, 2, 3, 4, 5, 6)
alpha   <- 0.05
init_seed <- 2021L      # fixed k-means init seed (shared by both methods)

n  <- n_per * k
true_clusters <- rep(seq_len(k), each = n_per)

# 3 cluster means on an equilateral triangle, scaled by delta
mu_for <- function(delta) rbind(c(delta / 2, 0),
                                c(0, sqrt(3) * delta / 2),
                                c(-delta / 2, 0))

# One replicate: returns a 1-row data.frame, or NULL on a hard error
# (e.g. k-means fails to return k distinct clusters). Warnings about
# non-exhaustive scanning are NOT errors -- the p-value is still valid.
one_rep <- function(r, delta, di) {
  mu <- mu_for(delta)
  set.seed(10000L + 1000L * di + r)                 # data seed
  X <- matrix(rnorm(n * q, sd = sig), n, q) + mu[true_clusters, ]

  fit <- tryCatch(
    suppressWarnings(suppressMessages(kmeans_inference_union(
      X, k = k, cluster_1 = 1, cluster_2 = 3,
      sig = sig, iter.max = iter.max, seed = init_seed,
      max_paths = max_paths, verbose = FALSE))),
    error = function(e) NULL)

  if (is.null(fit)) return(NULL)
  data.frame(
    delta      = delta,
    rep        = r,
    p_union    = fit$pval_union,
    p_path     = fit$pval_path,
    p_naive    = fit$p_naive,
    n_paths    = length(fit$paths),
    exhaustive = fit$exhaustive)
}

t0 <- Sys.time()
rows <- list()
for (di in seq_along(deltas)) {
  delta <- deltas[di]
  reps  <- mclapply(seq_len(nreps), one_rep, delta = delta, di = di,
                    mc.cores = ncores, mc.preschedule = FALSE)
  reps  <- Filter(Negate(is.null), reps)
  rows  <- c(rows, reps)
  cat(sprintf("delta=%.1f done (%d/%d valid) at %.0fs\n",
              delta, length(reps), nreps,
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

res <- do.call(rbind, rows)
saveRDS(res, "sims/results/power_type1_raw.rds")
write.csv(res, "sims/results/power_type1_raw.csv", row.names = FALSE)

# ---- summary: rejection rate at alpha by delta and method --------------------
agg <- aggregate(cbind(p_union, p_path, p_naive) ~ delta, data = res,
                 FUN = function(p) mean(p <= alpha))
agg_n <- aggregate(rep ~ delta, data = res, FUN = length)
names(agg_n)[2] <- "n_valid"
agg_paths <- aggregate(cbind(n_paths, exhaustive) ~ delta, data = res, FUN = mean)
summary_tab <- merge(merge(agg, agg_n), agg_paths)
names(summary_tab) <- c("delta", "rej_union", "rej_path", "rej_naive",
                        "n_valid", "mean_n_paths", "frac_exhaustive")

cat("\n==== Rejection rate at alpha =", alpha, "====\n")
print(summary_tab, row.names = FALSE, digits = 3)

# KS uniformity test under the global null (delta = 0) -------------------------
null0 <- res[res$delta == 0, ]
if (nrow(null0) > 0) {
  ks_u <- suppressWarnings(ks.test(null0$p_union, "punif"))
  ks_p <- suppressWarnings(ks.test(null0$p_path,  "punif"))
  cat(sprintf("\nGlobal-null (delta=0) uniformity KS p-value: union=%.3f  path=%.3f\n",
              ks_u$p.value, ks_p$p.value))
}

# power improvement check ------------------------------------------------------
pos <- summary_tab[summary_tab$delta > 0, ]
cat(sprintf("\nUnion >= path rejection at every delta>0: %s\n",
            all(pos$rej_union + 1e-12 >= pos$rej_path)))
cat(sprintf("Mean power gain (union - path) over delta>0: %+.3f\n",
            mean(pos$rej_union - pos$rej_path)))

saveRDS(summary_tab, "sims/results/power_type1_summary.rds")
cat("\nWrote sims/results/power_type1_{raw.csv,raw.rds,summary.rds}\n")
