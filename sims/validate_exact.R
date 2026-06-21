#!/usr/bin/env Rscript
# Calibration of the EXACT R-fiber unknown-variance union/path test.
# Global null (H0' holds), K=3, q in {2,10}, n=60. Should be ~Uniform.
suppressMessages(library(parallel))
source("sims/kmeans_union_unknownvar_exact.R")
nreps <- 200L; n_per <- 20L; k <- 3L; alpha <- 0.05
ncores <- if (length(commandArgs(TRUE))) as.integer(commandArgs(TRUE)[1]) else 4L
one <- function(rep, q) {
  set.seed(80000L + 1000L * q + rep)
  X <- matrix(rnorm(n_per * k * q), n_per * k, q)
  r <- tryCatch(kmeans_union_unknownvar_exact(X, k, 1, 3), error = function(e) NULL)
  if (is.null(r)) return(NULL)
  data.frame(q = q, p_union = r$p_union, p_path = r$p_path, n_int = r$n_intervals)
}
ks  <- function(p) { p <- p[is.finite(p)]; if (length(p) >= 2) suppressWarnings(ks.test(p, "punif")$p.value) else NA }
rej <- function(p) mean(p <= alpha, na.rm = TRUE)
cat("EXACT R-fiber unknown-variance UNION calibration, K=3 global null n=60\n")
for (q in c(2L, 10L)) {
  rr <- do.call(rbind, Filter(Negate(is.null),
        mclapply(seq_len(nreps), one, q = q, mc.cores = ncores, mc.preschedule = FALSE)))
  cat(sprintf("q=%-2d nval=%d mean#int=%.2f multi=%.2f | UNION rej=%.3f KS=%.3f | PATH rej=%.3f KS=%.3f\n",
    q, nrow(rr), mean(rr$n_int), mean(rr$n_int > 1),
    rej(rr$p_union), ks(rr$p_union), rej(rr$p_path), ks(rr$p_path)))
}
cat("DONE\n")
