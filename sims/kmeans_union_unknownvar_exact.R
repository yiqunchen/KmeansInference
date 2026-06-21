# ===========================================================================
# EXACT unknown-variance union/path p-value for k-means (any K), via the
# R-fiber arc-sweep. No Monte Carlo. Requires H0' (pointwise within-tested-
# cluster null). Replaces Yun-Barber's importance sampling for K>2.
#   - sweep theta in (0, pi/2); at each landing run deterministic Lloyd, get the
#     trajectory, compute its EXACT arc (Module 8), and if its output preserves
#     the tested clusters A,B add the arc to the union; jump past the arc.
#   - map arcs to r = d tan^2(theta); apply truncated-F_{q,(m-2)q}.
# ===========================================================================
suppressMessages(library(KmeansInference))
source("sims/unknown_var_union.R")        # fun_P0, fun_P1
source("sims/rfiber_exact.R")             # modular exact solver (+ unit tests)
.lloyd_  <- KmeansInference:::.kmeans_estimation_fixed_init
.presv_  <- KmeansInference:::.kmeans_preserves_observed_clusters

.build_fiber <- function(X, cl, k1, k2) {
  P0X <- fun_P0(cl, k1, k2) %*% X; P1X <- fun_P1(cl, k1, k2) %*% X; P2X <- X - P0X - P1X
  a0 <- sqrt(sum(P0X^2)); a1 <- sqrt(sum(P1X^2))
  if (a0 <= 0 || a1 <= 0) stop("degenerate P0X/P1X")
  Tt <- a0^2 + a1^2; d <- sum(cl == k1) + sum(cl == k2) - 2
  list(cmat = P2X, amat = sqrt(Tt) * (P0X / a0), bmat = sqrt(Tt) * (P1X / a1),
       d = d, R_obs = d * a0^2 / a1^2)
}

kmeans_union_unknownvar_exact <- function(X, k, cluster_1, cluster_2, seed = 2021,
                                          iter.max = 10, eps = 1e-7, max_steps = 5000L) {
  est <- kmeans_estimation(X, k, iter.max, seed, verbose = FALSE)
  cl <- est$final_cluster
  if (length(unique(cl)) < k) stop("k-means did not return k distinct clusters")
  init <- est$random_init_obs; k1 <- cluster_1; k2 <- cluster_2
  fib <- .build_fiber(X, cl, k1, k2); d <- fib$d; q <- ncol(X); df1 <- q; df2 <- d * q
  fl  <- list(cmat = fib$cmat, amat = fib$amat, bmat = fib$bmat)
  xth <- function(th) fib$cmat + sin(th) * fib$amat + cos(th) * fib$bmat

  union_arcs <- list(); th <- eps; hi <- pi / 2 - eps; steps <- 0L
  while (th < hi && steps < max_steps) {
    steps <- steps + 1L
    fr <- tryCatch(.lloyd_(xth(th), k, init, iter.max, tol_eps = 1e-6, verbose = FALSE),
                   error = function(e) NULL)
    if (is.null(fr)) { th <- th + 1e-4; next }            # degenerate slice; nudge
    traj <- fr$cluster[seq_len(fr$iter)]
    arc  <- path_arc(fl, traj, init, k)
    comp <- if (nrow(arc)) arc[arc[, 1] <= th + 1e-7 & th <= arc[, 2] + 1e-7, , drop = FALSE]
            else matrix(numeric(0), ncol = 2)
    if (nrow(comp) == 0) { th <- th + 1e-4; next }         # numerical fallback
    comp <- comp[1, ]
    if (.presv_(cl, fr$final_cluster, k1, k2)) union_arcs[[length(union_arcs) + 1]] <- comp
    th <- comp[2] + eps                                    # jump past this arc
  }
  if (length(union_arcs) == 0) stop("empty union (theta_obs should be preserved)")
  arcs  <- .merge_intervals(do.call(rbind, union_arcs))
  ivl_r <- cbind(r_of_theta(arcs[, 1], d), r_of_theta(arcs[, 2], d))
  R_obs <- fib$R_obs
  contains <- ivl_r[, 1] <= R_obs & R_obs <= ivl_r[, 2]
  ivl_path <- if (any(contains)) ivl_r[contains, , drop = FALSE]
              else ivl_r[which.min(abs(rowMeans(ivl_r) - R_obs)), , drop = FALSE]
  list(p_union = trunc_F_p(R_obs, ivl_r, df1, df2),
       p_path  = trunc_F_p(R_obs, ivl_path, df1, df2),
       R_obs = R_obs, n_intervals = nrow(ivl_r), n_steps = steps)
}
