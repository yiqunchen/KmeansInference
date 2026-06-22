# Plug-in variance (Chen-Witten) union/path p-values, computed CHEAPLY:
# the known-variance union test uses a sigma-INDEPENDENT truncation set on the
# phi-line. We get that set with a cheap phi-grid scan (like the R-fiber numerical
# version), then evaluate the chi-ratio p-value for any sigma (true, MED, Sample).
# Verified to reproduce KmeansInference::kmeans_inference_union for sigma_true.
suppressMessages(library(KmeansInference))
.perturb <- KmeansInference:::.kmeans_perturb_data
.lloyd   <- KmeansInference:::.kmeans_estimation_fixed_init
.presv   <- KmeansInference:::.kmeans_preserves_observed_clusters
.ivpval  <- KmeansInference:::.kmeans_interval_pvalue
.med     <- KmeansInference:::.kmeans_estimate_MED

# returns list with p-values (union & path) for sigma in c(known, MED, Sample)
kmeans_plugin_unknownvar <- function(X, k, c1, c2, sig_true = 1, seed = 2021,
                                     iter.max = 10, n_phi = 400, tol = 1e-8) {
  est <- kmeans_estimation(X, k, iter.max, seed, verbose = FALSE); cl <- est$final_cluster
  if (length(unique(cl)) < k) stop("k-means did not return k distinct clusters")
  init <- est$random_init_obs; q <- ncol(X); n <- nrow(X)
  n1 <- sum(cl == c1); n2 <- sum(cl == c2); sqnu <- 1/n1 + 1/n2; vn <- sqrt(sqnu)
  v_vec <- numeric(n); v_vec[cl == c1] <- 1/n1; v_vec[cl == c2] <- -1/n2
  diffm <- colMeans(X[cl == c1, , drop = FALSE]) - colMeans(X[cl == c2, , drop = FALSE])
  phi_obs <- sqrt(sum(diffm^2)); dir_XTv <- diffm / phi_obs

  s_med  <- .med(X)
  s_samp <- sqrt(sum(scale(X, center = TRUE, scale = FALSE)^2) / (length(X) - q))
  preserved <- function(phi) {
    Xp <- .perturb(X, phi, phi_obs, dir_XTv, v_vec, vn)
    fr <- tryCatch(.lloyd(Xp, k, init, iter.max, tol_eps = 1e-6, verbose = FALSE), error = function(e) NULL)
    if (is.null(fr)) FALSE else isTRUE(.presv(cl, fr$final_cluster, c1, c2))
  }
  refine <- function(t_false, t_true) { for (i in 1:50) { m <- (t_false + t_true)/2
    if (preserved(m)) t_true <- m else t_false <- m; if (abs(t_true - t_false) < tol) break }; (t_false + t_true)/2 }
  # scan phi up to where the chi tail is negligible for the LARGEST sigma used
  phi_max <- phi_obs + 8 * vn * max(sig_true, s_med, s_samp)
  ph <- seq(1e-7, phi_max, length.out = n_phi); pr <- vapply(ph, preserved, logical(1))
  if (!any(pr)) stop("empty preservation set")
  runs <- which(diff(c(FALSE, pr, FALSE)) != 0)
  st <- runs[seq(1, length(runs), 2)]; en <- runs[seq(2, length(runs), 2)] - 1L
  ivl <- t(mapply(function(s, e) c(
    if (s == 1L) ph[1] else refine(ph[s-1L], ph[s]),
    if (e == n_phi) ph[n_phi] else refine(ph[e+1L], ph[e])), st, en))
  if (is.null(dim(ivl))) ivl <- matrix(ivl, ncol = 2)
  Su <- intervals::reduce(intervals::Intervals_full(ivl), check_valid = FALSE)
  cont <- ivl[, 1] <= phi_obs & phi_obs <= ivl[, 2]
  Sp <- intervals::reduce(intervals::Intervals_full(
    if (any(cont)) ivl[cont, , drop = FALSE]
    else ivl[which.min(abs(rowMeans(ivl) - phi_obs)), , drop = FALSE]), check_valid = FALSE)
  pv <- function(S, sig) tryCatch(.ivpval(S, phi_obs, sqnu * sig^2, q), error = function(e) NA_real_)
  list(p_union_known = pv(Su, sig_true), p_path_known = pv(Sp, sig_true),
       p_union_med = pv(Su, s_med),  p_path_med  = pv(Sp, s_med),
       p_union_samp = pv(Su, s_samp), p_path_samp = pv(Sp, s_samp),
       sig_med = s_med, sig_samp = s_samp, n_int = nrow(ivl))
}
