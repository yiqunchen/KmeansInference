#!/usr/bin/env Rscript
# Validate the R-fiber unknown-variance UNION test: global null (H0' holds),
# K=3, q in {2,10}, n=60. The union p-value should be ~Uniform (rej ~ alpha,
# KS not tiny). Path is the single-component anchor (must also calibrate).
suppressMessages(library(parallel))
source("sims/kmeans_union_unknownvar.R")
nreps <- 150L; n_per <- 20L; k <- 3L; alpha <- 0.05
ncores <- if (length(commandArgs(TRUE))) as.integer(commandArgs(TRUE)[1]) else 4L

one <- function(rep, q) {
  set.seed(80000L + 1000L * q + rep)
  X <- matrix(rnorm(n_per * k * q), n_per * k, q)          # global null mu=0
  r <- tryCatch(kmeans_union_unknownvar(X, k = k, cluster_1 = 1, cluster_2 = 3,
                                        n_theta = 250, tol = 1e-7),
                error = function(e) NULL)
  if (is.null(r)) return(NULL)
  data.frame(q = q, p_union = r$p_union, p_path = r$p_path, n_int = r$n_intervals)
}
ks  <- function(p) { p <- p[is.finite(p)]; if (length(p) >= 2) suppressWarnings(ks.test(p, "punif")$p.value) else NA }
rej <- function(p) mean(p <= alpha, na.rm = TRUE)
cat("R-FIBER unknown-variance UNION validation, K=3 global null n=60\n")
for (q in c(2L, 10L)) {
  rr <- do.call(rbind, Filter(Negate(is.null),
        mclapply(seq_len(nreps), one, q = q, mc.cores = ncores, mc.preschedule = FALSE)))
  cat(sprintf("q=%-2d nval=%d mean#int=%.2f frac_multi=%.2f | UNION rej=%.3f KS=%.3f | PATH rej=%.3f KS=%.3f\n",
    q, nrow(rr), mean(rr$n_int), mean(rr$n_int > 1),
    rej(rr$p_union), ks(rr$p_union), rej(rr$p_path), ks(rr$p_path)))
}
cat("DONE\n")
