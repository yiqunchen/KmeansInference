# ===========================================================================
# VALID more-powerful UNION test under UNKNOWN variance (Yun-Barber F-pivot),
# computed on the R-FIBER (not by remapping the phi-fiber).
#
# Construction (per the proposition):
#   Fix the observed tested clusters A=hat C_a, B=hat C_b and the projections
#   P0 (between, rank 1), P1 (within-tested-cluster, rank d=m-2), P2=I-P0-P1.
#   Condition on Z = (T, U0=P0X/||P0X||, U1=P1X/||P1X||, P2X), T=||P0X||^2+||P1X||^2.
#   The R-fiber is x_z(theta) = P2X + sqrt(T)(sin(theta) U0 + cos(theta) U1),
#   theta in (0, pi/2), with R = d tan^2(theta) the only free coordinate; under
#   H0' (pointwise within-tested-cluster null), R ~ F_{q, dq} before selection.
#   Selection event E_{A,B} = { deterministic Lloyd run on x_z outputs A and B }.
#   S'_union(Z) = { theta : x_z(theta) in E_{A,B} }  ->  R-set via r = d tan^2 theta.
#   p = P(R >= R_obs, R in S'_union) / P(R in S'_union),  R ~ F_{q, dq}.
#
# This file finds S'_union NUMERICALLY (fine theta sweep + bisection-refined
# boundaries, exact to ~tol). The exact route is quartic-in-tan(theta/2); the
# numerical version validates the statistical claim first.
# ===========================================================================
suppressMessages(library(KmeansInference))
source("sims/unknown_var_union.R")               # fun_P0, fun_P1, fun_ts
.lloyd  <- KmeansInference:::.kmeans_estimation_fixed_init
.presrv <- KmeansInference:::.kmeans_preserves_observed_clusters

# truncated-F upper p-value over a set of r-intervals (rows [lo,hi]).
.trunc_F_p <- function(R_obs, ivl, df1, df2) {
  if (is.null(ivl) || nrow(ivl) == 0) return(NA_real_)
  pint <- function(a, b) pf(b, df1, df2) - pf(a, df1, df2)
  den <- sum(apply(ivl, 1, function(z) pint(z[1], z[2]))); if (den <= 0) return(NA_real_)
  num <- sum(apply(ivl, 1, function(z) { a <- max(z[1], R_obs); if (z[2] <= a) 0 else pint(a, z[2]) }))
  max(0, min(1, num / den))
}

#' Unknown-variance union (and path) p-values on the R-fiber.
#' @return list(p_union, p_path, R_obs, n_intervals, intervals_r, theta_obs)
kmeans_union_unknownvar <- function(X, k, cluster_1, cluster_2, seed = 2021,
                                    iter.max = 10, n_theta = 300, tol = 1e-9) {
  est <- kmeans_estimation(X, k, iter.max, seed, verbose = FALSE)
  cl  <- est$final_cluster
  if (length(unique(cl)) < k) stop("k-means did not return k distinct clusters")
  init <- est$random_init_obs
  k1 <- cluster_1; k2 <- cluster_2
  m <- sum(cl == k1) + sum(cl == k2); d <- m - 2
  if (d < 1) stop("need m > 2")
  P0X <- fun_P0(cl, k1, k2) %*% X; P1X <- fun_P1(cl, k1, k2) %*% X; P2X <- X - P0X - P1X
  a0 <- sqrt(sum(P0X^2)); a1 <- sqrt(sum(P1X^2))
  if (a0 <= 0 || a1 <= 0) stop("degenerate P0X/P1X")
  Tt <- a0^2 + a1^2; q <- ncol(X); df1 <- q; df2 <- d * q
  A2 <- sqrt(Tt) * (P0X / a0); B2 <- sqrt(Tt) * (P1X / a1)   # xz(th)=P2X+sin*A2+cos*B2
  R_obs <- d * a0^2 / a1^2
  theta_obs <- atan(sqrt(R_obs / d))

  preserved <- function(theta) {
    Xz <- P2X + sin(theta) * A2 + cos(theta) * B2
    fr <- tryCatch(.lloyd(Xz, k, init, iter.max, tol_eps = 1e-6, verbose = FALSE),
                   error = function(e) NULL)
    if (is.null(fr)) return(FALSE)
    isTRUE(.presrv(cl, fr$final_cluster, k1, k2))
  }
  # find the FALSE/TRUE boundary in [t_false, t_true] by bisection
  refine <- function(t_false, t_true) {
    for (i in 1:60) { mid <- (t_false + t_true) / 2
      if (preserved(mid)) t_true <- mid else t_false <- mid
      if (abs(t_true - t_false) < tol) break }
    (t_false + t_true) / 2
  }

  th  <- seq(1e-7, pi/2 - 1e-7, length.out = n_theta)
  pr  <- vapply(th, preserved, logical(1))
  if (!any(pr)) stop("empty preservation set (theta_obs should be preserved)")
  # extract maximal TRUE runs, refine the two boundaries of each
  runs <- which(diff(c(FALSE, pr, FALSE)) != 0)
  starts <- runs[seq(1, length(runs), 2)]; ends <- runs[seq(2, length(runs), 2)] - 1L
  ivl_theta <- t(mapply(function(s, e) {
    lo <- if (s == 1L) th[1] else refine(th[s - 1L], th[s])     # FALSE@s-1 -> TRUE@s
    hi <- if (e == n_theta) th[n_theta] else refine(th[e + 1L], th[e]) # FALSE@e+1 -> TRUE@e
    c(lo, hi)
  }, starts, ends))
  if (is.null(dim(ivl_theta))) ivl_theta <- matrix(ivl_theta, ncol = 2)
  ivl_r <- cbind(d * tan(ivl_theta[, 1])^2, d * tan(ivl_theta[, 2])^2)

  # path = the interval containing R_obs (connected component)
  contains <- ivl_r[, 1] <= R_obs & R_obs <= ivl_r[, 2]
  ivl_path <- if (any(contains)) ivl_r[contains, , drop = FALSE]
              else ivl_r[which.min(abs(rowMeans(ivl_r) - R_obs)), , drop = FALSE]

  list(p_union = .trunc_F_p(R_obs, ivl_r,   df1, df2),
       p_path  = .trunc_F_p(R_obs, ivl_path, df1, df2),
       R_obs = R_obs, n_intervals = nrow(ivl_r), intervals_r = ivl_r,
       theta_obs = theta_obs)
}
