#!/usr/bin/env Rscript
# ===========================================================================
# exp_datathin.R -- compare DATA THINNING (Neufeld-Gao-Popp-Witten) against the
# selective k-means tests on the SAME data, per-rep matched, for three metrics:
#   (1) Type I error      -- global null (single blob), test a split pair.
#   (2) detection prob     -- does the clustering recover the true pair?
#   (3) conditional power  -- rejection | the pair was recovered.
#
# Both sigma-regimes, matched to the fully-conditional case:
#   KNOWN sigma   : selective UNION         vs  thinning + TRUE sigma^2 (chi^2 on test fold)
#   UNKNOWN sigma : selective STUDENTIZED-F vs  thinning + sigma^2_MED plug-in, F-test on fold
#   naive z/chi^2 = invalid reference.
#
# Data thinning (isotropic N(mu, sigma^2 I), epsilon=0.5): datathin entrywise
#   X = X1 + X2, X1,X2 indep ~ N(.5 mu, .5 sigma^2 I) AT THE VARIANCE USED TO SPLIT.
#   Cluster X1, TEST on X2 (X2 _|_ X1 => plain chi^2/F mean test is valid, no
#   conditioning). KEY ASYMMETRY: thinning needs sigma^2 to SPLIT (true=oracle,
#   estimate=plug-in); the selective studentized test needs NO sigma at all.
#   The tutorial's gamma route (family='normal-variance') is the ORTHOGONAL
#   problem (mean known, infer variance); the gamma/chi^2 reappears here only as
#   the F-test denominator (pooled within-cluster SS), which makes the TEST sigma-free.
# Saves per-rep rows -> sims/results/datathin_compare.rds  (plot_datathin.R draws).
# ===========================================================================
suppressMessages({ library(KmeansInference); library(datathin); library(parallel) })
source("sims/unknown_var_union.R")
source("sims/kmeans_union_unknownvar.R")
MED <- KmeansInference:::.kmeans_estimate_MED              # robust SD estimator (Chen-Witten)

.envi <- function(k, d) { v <- Sys.getenv(k); if (nzchar(v)) as.integer(v) else d }
K <- 3L; Q <- 2L; SIG <- 1; NPER <- 10L; EPS <- 0.5
INIT_SEED <- 2021L; ITERMAX <- 10L; MAXPATHS <- 100L; STP <- 1e-6; NTHETA <- 250L
NREP_NULL <- .envi("DT_NULL", 400L); NREP_SIG <- .envi("DT_SIG", 220L)
DELTAS <- if (nzchar(Sys.getenv("DT_DELTAS"))) as.numeric(strsplit(Sys.getenv("DT_DELTAS"), ",")[[1]]) else c(4, 5, 6, 7)
NCORES <- .envi("DT_CORES", max(1L, detectCores() - 2L))

tri_means <- function(delta) rbind(c(delta/2, 0), c(0, sqrt(3)*delta/2), c(-delta/2, 0))

## recover a target TRUE pair from a clustering (>=75% purity, distinct) ------
recover_pair <- function(cl, true, target = c(1L, 2L), pur = 0.75) {
  find <- function(tc) {
    best <- NA_integer_; bestpur <- 0
    for (a in unique(cl)) {
      tt <- true[cl == a]; if (!length(tt)) next
      tb <- table(tt); dom <- as.integer(names(which.max(tb))); pp <- max(tb)/sum(tb)
      if (dom == tc && pp >= pur && pp > bestpur) { best <- a; bestpur <- pp }
    }
    best
  }
  a <- find(target[1]); b <- find(target[2])
  if (is.na(a) || is.na(b) || a == b) return(NULL)
  c(a, b)
}

## thinning: cluster X1, test the cluster-mean difference on X2 ----------------
# split_var = variance handed to datathin (true sigma^2 = oracle, or sigma^2_MED).
# test = 'chi2' (oracle known sigma, tau^2 = eps*split_var) or 'F' (sigma-free).
thin_once <- function(X, true, is_null, split_var, test, sigma_for_chi2 = NA) {
  th <- tryCatch(datathin(X, family = "normal", arg = split_var, epsilon = c(EPS, 1 - EPS)),
                 error = function(e) NULL)
  if (is.null(th)) return(list(detected = FALSE, p = NA_real_))
  X1 <- th[, , 1]; X2 <- th[, , 2]
  f1 <- tryCatch(kmeans_estimation(X1, K, ITERMAX, INIT_SEED, verbose = FALSE), error = function(e) NULL)
  if (is.null(f1)) return(list(detected = FALSE, p = NA_real_))
  cl1 <- f1$final_cluster
  pr <- if (is_null) c(1L, 2L) else recover_pair(cl1, true)
  if (is.null(pr)) return(list(detected = FALSE, p = NA_real_))
  Ca <- which(cl1 == pr[1]); Cb <- which(cl1 == pr[2]); na <- length(Ca); nb <- length(Cb)
  d   <- colMeans(X2[Ca, , drop = FALSE]) - colMeans(X2[Cb, , drop = FALSE])
  ssd <- sum(d^2) / (1/na + 1/nb)                          # ~ tau^2 chi^2_q under H0
  if (test == "chi2") {
    tau2 <- EPS * sigma_for_chi2^2
    p <- stats::pchisq(ssd / tau2, Q, lower.tail = FALSE)
  } else {                                                 # self-calibrating F (sigma-free)
    ssw <- 0; dfw <- 0L
    for (a in unique(cl1)) { idx <- which(cl1 == a)
      if (length(idx) >= 2) { cm <- colMeans(X2[idx, , drop = FALSE])
        ssw <- ssw + sum(sweep(X2[idx, , drop = FALSE], 2, cm)^2); dfw <- dfw + (length(idx) - 1L) * Q } }
    p <- if (dfw < 1 || ssw <= 0) NA_real_ else stats::pf((ssd/Q)/(ssw/dfw), Q, dfw, lower.tail = FALSE)
  }
  list(detected = TRUE, p = p)
}

## one replicate --------------------------------------------------------------
one_rep <- function(r, delta, is_null, seed0) {
  n <- NPER * K; true <- rep(seq_len(K), each = NPER)
  mu <- if (is_null) matrix(0, K, Q) else tri_means(delta)
  set.seed(seed0 + r)
  X <- matrix(rnorm(n * Q), n, Q) * SIG + mu[true, , drop = FALSE]
  sig_med  <- tryCatch(MED(X), error = function(e) SIG)                                  # robust (median) SD
  sig_samp <- sqrt(sum(scale(X, center = TRUE, scale = FALSE)^2) / (length(X) - Q))      # total sample SD

  out <- data.frame(rep = r, delta = delta, is_null = is_null,
                    sel_detected = FALSE, thinO_detected = FALSE, thinM_detected = FALSE, thinS_detected = FALSE,
                    sig_med = sig_med, sig_samp = sig_samp,
                    p_naive = NA_real_, p_path = NA_real_, p_union = NA_real_, p_studF = NA_real_,
                    p_thin_oracle = NA_real_, p_thin_med = NA_real_, p_thin_samp = NA_real_)

  ## SELECTIVE: cluster the FULL data
  full <- tryCatch(kmeans_estimation(X, K, ITERMAX, INIT_SEED, verbose = FALSE), error = function(e) NULL)
  if (!is.null(full)) {
    pr <- if (is_null) c(1L, 2L) else recover_pair(full$final_cluster, true)
    if (!is.null(pr)) {
      out$sel_detected <- TRUE
      f <- tryCatch(suppressWarnings(suppressMessages(kmeans_inference_union(
             X, k = K, cluster_1 = pr[1], cluster_2 = pr[2], iso = TRUE, sig = SIG,
             iter.max = ITERMAX, seed = INIT_SEED, max_paths = MAXPATHS,
             scan_tail_prob = STP, verbose = FALSE))), error = function(e) NULL)
      if (!is.null(f)) { out$p_naive <- f$p_naive; out$p_path <- f$pval_path; out$p_union <- f$pval_union }
      rf <- tryCatch(kmeans_union_unknownvar(X, k = K, cluster_1 = pr[1], cluster_2 = pr[2],
              seed = INIT_SEED, n_theta = NTHETA, tol = 1e-7), error = function(e) NULL)
      if (!is.null(rf)) out$p_studF <- rf$p_union
    }
  }

  ## THINNING -- known sigma (true sigma^2 split, chi^2 test)
  to <- thin_once(X, true, is_null, split_var = SIG^2,      test = "chi2", sigma_for_chi2 = SIG)
  out$thinO_detected <- to$detected; out$p_thin_oracle <- to$p
  ## THINNING -- unknown sigma, sigma^2_MED split, self-calibrating F-test
  tm <- thin_once(X, true, is_null, split_var = sig_med^2,  test = "F")
  out$thinM_detected <- tm$detected; out$p_thin_med <- tm$p
  ## THINNING -- unknown sigma, sample-variance split, self-calibrating F-test
  ts <- thin_once(X, true, is_null, split_var = sig_samp^2, test = "F")
  out$thinS_detected <- ts$detected; out$p_thin_samp <- ts$p
  out
}

run_cell <- function(delta, is_null, seed0, nrep, label) {
  cat(sprintf("[%s] delta=%g null=%s nrep=%d ...\n", label, delta, is_null, nrep))
  rows <- mclapply(seq_len(nrep), function(r) tryCatch(one_rep(r, delta, is_null, seed0),
                   error = function(e) NULL), mc.cores = NCORES, mc.preschedule = FALSE)
  do.call(rbind, Filter(Negate(is.null), rows))
}

t0 <- Sys.time()
cells <- list()
cells[["null"]] <- run_cell(0, TRUE, 30000L, NREP_NULL, "TypeI null")
for (i in seq_along(DELTAS)) { d <- DELTAS[i]
  cells[[paste0("d", d)]] <- run_cell(d, FALSE, 40000L + i * 1000L, NREP_SIG, paste0("signal d=", d)) }
all <- do.call(rbind, cells)
saveRDS(all, "sims/results/datathin_compare.rds")

rej <- function(p) mean(p <= 0.05, na.rm = TRUE)
cat("\n==== Type I (null cell) ====\n")
nu <- cells[["null"]]
cat(sprintf("naive=%.3f path=%.3f union=%.3f studF=%.3f | thin_oracle=%.3f thin_MED=%.3f thin_samp=%.3f\n",
            rej(nu$p_naive), rej(nu$p_path), rej(nu$p_union), rej(nu$p_studF),
            rej(nu$p_thin_oracle), rej(nu$p_thin_med), rej(nu$p_thin_samp)))
cat(sprintf("   (null sig estimates: MED=%.2f  sample=%.2f  true=%.2f)\n", median(nu$sig_med), median(nu$sig_samp), SIG))
cat("\n==== detection prob (full vs thinned) & conditional power ====\n")
for (d in DELTAS) { cc <- cells[[paste0("d", d)]]
  cp <- function(p, det) { m <- det & is.finite(p); if (!sum(m)) NA else mean(p[m] <= 0.05) }
  cat(sprintf("d=%g | sigMED=%.2f sigSamp=%.2f | detect sel=%.2f thinO=%.2f thinM=%.2f thinS=%.2f | cpow union=%.2f studF=%.2f thinO=%.2f thinM=%.2f thinS=%.2f\n",
      d, median(cc$sig_med), median(cc$sig_samp), mean(cc$sel_detected),
      mean(cc$thinO_detected), mean(cc$thinM_detected), mean(cc$thinS_detected),
      cp(cc$p_union, cc$sel_detected), cp(cc$p_studF, cc$sel_detected),
      cp(cc$p_thin_oracle, cc$thinO_detected), cp(cc$p_thin_med, cc$thinM_detected), cp(cc$p_thin_samp, cc$thinS_detected)))
}
cat(sprintf("\nDone in %.1f min. Wrote sims/results/datathin_compare.rds (%d rows)\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), nrow(all)))
