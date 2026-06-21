#!/usr/bin/env Rscript
# ===========================================================================
# Generalized sweep harness for the union vs. path selective k-means test.
#
# Goal
# ----
# Run ONE parameter cell (run_cell) and return a one-row per-cell summary that,
# across many cells, establishes that the union test
#   (a) controls selective Type I error at level alpha (uniform null p-values),
#   (b) is more powerful than the original path-conditioned test, with the
#       naive z/chisq test shown as an invalid (over-rejecting) reference.
#
# It is built on sims/power_type1_union_vs_path.R: same parallel
# mclapply-over-reps structure, and the SAME "only skip on hard error, never on
# warning" logic (a warning about non-exhaustive scanning does NOT invalidate
# the p-value; only a hard error -- e.g. k-means failing to return k distinct
# clusters -- drops a rep).
#
# Dependency-light: base R + KmeansInference (+ parallel, which ships with R).
#
# ---------------------------------------------------------------------------
# CONFIG INTERFACE: a "cell" is a named list cfg with these fields (defaults
# in parentheses; .default_cfg() supplies them, .fill_cfg() merges yours in).
#
#   --- data dimensions / cluster geometry ---
#   n_per        (8)      observations per TRUE cluster, used when `balance`
#                         is NULL (the balanced case: n = n_per * k).
#   n            (NULL)   total n; used WITH `balance` for unbalanced designs.
#   balance      (NULL)   length-k numeric weights; cluster sizes are
#                         round(n * balance / sum(balance)). Overrides n_per.
#   k            (3)      number of clusters for k-means.
#   q            (2)      data dimension.
#   geometry     ('triangle')  true mean configuration:
#                  'triangle'    : k means on a regular simplex / for k=3 an
#                                  equilateral triangle, scaled by delta. The
#                                  canonical separated-clusters design.
#                  'single-blob' : all means at the origin (a pure global null;
#                                  delta is ignored). Use for Type-I cells.
#                  'elongated'   : 'triangle' geometry but the means are
#                                  stretched along coordinate 1 by a factor
#                                  `elong` (default 3), giving anisotropic
#                                  separation (stresses the covariance path).
#   delta        (0)      separation scale of the cluster means. delta = 0 is
#                         the global null for 'triangle'/'elongated'.
#   elong        (3)      stretch factor for geometry = 'elongated'.
#
#   --- noise / covariance model ---
#   sig          (NULL)   isotropic noise SD used to SIMULATE the data. If a
#                         number, that SD is used for generation. The value
#                         passed to the TEST is `sig_test` (see below).
#   sig_test     (sig)    isotropic SD handed to the iso test. NULL => the test
#                         estimates it with the robust median estimator. Set
#                         this != sig to probe robustness to a misspecified SD.
#   Sigma        (NULL)   q x q SPD covariance used to SIMULATE the data. If
#                         non-NULL the general-covariance path is used:
#                         iso = FALSE and the test receives SigInv_test.
#   SigInv_test  (solve(Sigma))  the inverse covariance ACTUALLY passed to the
#                         test. Default is the correct inverse; set to a wrong
#                         matrix to probe misspecification. Ignored if Sigma is
#                         NULL.  (When Sigma is set, sig/sig_test are ignored.)
#
#   --- the tested contrast ---
#   pair         (c(1,3))   c(cluster_1, cluster_2): the two OBSERVED k-means
#                           cluster labels whose mean difference is tested.
#   tested_pair_type ('different')  bookkeeping label echoed into the summary,
#                           describing whether the tested pair corresponds to a
#                           genuine mean difference ('different') or to a split
#                           of a single true blob ('null-split'). It does NOT
#                           change the data; choose geometry/delta to realize it
#                           (e.g. 'single-blob' + 'null-split', or
#                           'triangle' delta>0 + 'different').
#
#   --- k-means / exploration controls (forwarded to the test) ---
#   iter.max       (10)
#   max_paths      (100)
#   scan_tail_prob (1e-6)
#
#   --- experiment controls ---
#   init_seed    (2021)   fixed k-means INIT seed, shared across reps & methods.
#   data_seed0   (10000)  base for per-rep DATA seeds (independent of init_seed).
#                         Rep r uses set.seed(data_seed0 + cell_seed_offset + r).
#   cell_seed_offset (0)  add a distinct offset per cell so different cells draw
#                         independent data even with the same data_seed0.
#   nreps        (30)
#   alpha        (0.05)
#   ncores       (1)      mclapply cores. KEEP AT 1 when run_cell is itself
#                         called from a parallel orchestrator (the default).
#   label        (NULL)   optional human-readable cell name, echoed in summary.
#
# A cell is treated as NULL (for Type-I / KS diagnostics) iff
#   geometry == 'single-blob'  OR  delta == 0  OR  tested_pair_type=='null-split'.
# For null cells the summary reports KS-uniformity p-values; for alt cells it
# reports rejection (power) and the union-minus-path power gain.
#
# ---------------------------------------------------------------------------
# OUTPUT of run_cell(cfg): a one-row data.frame with the echoed cfg plus:
#   n_valid, rej_union, rej_path, rej_naive (mean p <= alpha),
#   ks_union, ks_path (KS-vs-uniform p; NA for non-null cells),
#   power_gain (mean(rej_union - rej_path); for null cells = type-I gap),
#   mean_n_paths, frac_exhaustive, median_runtime_s, is_null, n_total.
#
# DRIVER: run_grid(cells, ...) maps run_cell over a list of cfgs and rbinds the
# summaries (its own ncores parallelizes ACROSS cells; each cell stays serial).
#
# Usage:
#   source("sims/sweep_harness.R")
#   run_cell(list(geometry="single-blob", nreps=30, n_per=8))
#   run_grid(list(cellA, cellB), grid_ncores = 4)
# Or as a script (runs a tiny 2-cell smoke grid):
#   Rscript sims/sweep_harness.R
# ===========================================================================

suppressMessages(library(KmeansInference))
suppressMessages(library(parallel))
# unknown-variance: (a) the valid PATH p-value via the phi-remap (post-processing
# of the fit), and (b) the VALID MORE-POWERFUL UNION computed on the R-fiber
# (Yun-Barber F-pivot in the pivot's own coordinate; the phi-remap union was
# anti-conservative). The R-fiber version re-explores along x'(r) ~4s/rep.
source("sims/unknown_var_union.R")
source("sims/kmeans_union_unknownvar.R")   # numerical R-fiber union/path

# ---- config defaults & merge ----------------------------------------------
.default_cfg <- function() list(
  n_per = 8L, n = NULL, balance = NULL,
  k = 3L, q = 2L,
  geometry = "triangle", delta = 0, elong = 3,
  sig = NULL, sig_test = NULL,
  Sigma = NULL, SigInv_test = NULL,
  pair = c(1L, 3L), tested_pair_type = "different",
  iter.max = 10L, max_paths = 100L, scan_tail_prob = 1e-6,
  init_seed = 2021L, data_seed0 = 10000L, cell_seed_offset = 0L,
  nreps = 30L, alpha = 0.05, ncores = 1L, label = NULL
)

.fill_cfg <- function(cfg) {
  d <- .default_cfg()
  for (nm in names(cfg)) d[[nm]] <- cfg[[nm]]
  # sig_test defaults to sig (the simulation SD) when not given explicitly
  if (is.null(cfg$sig_test) && !("sig_test" %in% names(cfg))) d$sig_test <- d$sig
  d
}

# ---- true cluster means by geometry ---------------------------------------
# Returns a k x q matrix of cluster means (rows = clusters).
.geometry_means <- function(geometry, k, q, delta, elong) {
  M <- matrix(0, nrow = k, ncol = q)
  if (geometry == "single-blob") return(M)   # all means at origin: global null

  if (geometry == "triangle" || geometry == "elongated") {
    if (k == 3L && q >= 2L) {
      base <- rbind(c(delta / 2, 0),
                    c(0, sqrt(3) * delta / 2),
                    c(-delta / 2, 0))
      M[, 1:2] <- base
    } else {
      # general k: place means on the vertices of a scaled regular simplex,
      # embedded in the first (k-1) coords (capped at q).
      d_eff <- min(k - 1L, q)
      V <- diag(1, nrow = k, ncol = d_eff)             # k x d_eff, one-hot
      V <- scale(V, center = TRUE, scale = FALSE)      # center the simplex
      # normalize so adjacent-vertex distance is delta
      if (k >= 2L) {
        pd <- sqrt(sum((V[1, ] - V[2, ])^2))
        if (pd > 0) V <- V * (delta / pd)
      }
      M[, seq_len(d_eff)] <- V
    }
    if (geometry == "elongated") M[, 1] <- M[, 1] * elong
    return(M)
  }
  stop(sprintf("unknown geometry '%s'", geometry))
}

# ---- true cluster sizes ----------------------------------------------------
.cluster_sizes <- function(cfg) {
  if (!is.null(cfg$balance)) {
    if (length(cfg$balance) != cfg$k) stop("balance must have length k")
    if (is.null(cfg$n)) stop("n must be supplied when balance is used")
    w  <- cfg$balance / sum(cfg$balance)
    sz <- round(cfg$n * w)
    # fix rounding so sizes sum to n
    diff <- cfg$n - sum(sz)
    if (diff != 0) sz[which.max(sz)] <- sz[which.max(sz)] + diff
    if (any(sz < 1)) stop("balance produced an empty cluster")
    return(as.integer(sz))
  }
  rep(as.integer(cfg$n_per), cfg$k)
}

.is_null_cell <- function(cfg) {
  cfg$geometry == "single-blob" ||
    isTRUE(cfg$delta == 0) ||
    identical(cfg$tested_pair_type, "null-split")
}

# ---- one replicate ---------------------------------------------------------
# Returns a 1-row data.frame, or NULL on a HARD error only. Warnings (e.g.
# non-exhaustive scanning) are suppressed and do NOT drop the rep.
.one_rep <- function(r, cfg, sizes, mu, use_genCov, SigChol) {
  n <- sum(sizes)
  true_clusters <- rep(seq_len(cfg$k), times = sizes)

  set.seed(cfg$data_seed0 + cfg$cell_seed_offset + r)   # per-rep DATA seed
  Z <- matrix(rnorm(n * cfg$q), n, cfg$q)
  if (use_genCov) {
    noise <- Z %*% SigChol                 # rows ~ N(0, Sigma)
  } else {
    sd_sim <- if (is.null(cfg$sig)) 1 else cfg$sig
    noise <- Z * sd_sim
  }
  X <- noise + mu[true_clusters, , drop = FALSE]

  t_start <- proc.time()[["elapsed"]]
  fit <- tryCatch(
    suppressWarnings(suppressMessages(
      if (use_genCov) {
        kmeans_inference_union(
          X, k = cfg$k, cluster_1 = cfg$pair[1], cluster_2 = cfg$pair[2],
          iso = FALSE, SigInv = cfg$SigInv_test,
          iter.max = cfg$iter.max, seed = cfg$init_seed,
          max_paths = cfg$max_paths, scan_tail_prob = cfg$scan_tail_prob,
          verbose = FALSE)
      } else {
        kmeans_inference_union(
          X, k = cfg$k, cluster_1 = cfg$pair[1], cluster_2 = cfg$pair[2],
          iso = TRUE, sig = cfg$sig_test,
          iter.max = cfg$iter.max, seed = cfg$init_seed,
          max_paths = cfg$max_paths, scan_tail_prob = cfg$scan_tail_prob,
          verbose = FALSE)
      }
    )),
    error = function(e) NULL)
  runtime <- proc.time()[["elapsed"]] - t_start

  if (is.null(fit)) return(NULL)

  # Detection (Chen & Witten eq. 24): are the two tested ESTIMATED clusters
  # genuine recovered true clusters? We call a tested estimated cluster
  # "recovered" if >= 75% of its members come from a single true cluster, and
  # the two tested clusters' dominant true clusters differ.
  est <- fit$final_cluster
  dom_pur <- function(lbl) {
    tt <- true_clusters[est == lbl]
    if (length(tt) == 0) return(c(dom = NA_real_, pur = 0))
    tb <- table(tt); c(dom = as.numeric(names(which.max(tb))), pur = max(tb) / sum(tb))
  }
  da <- dom_pur(cfg$pair[1]); db <- dom_pur(cfg$pair[2])
  detected <- is.finite(da["dom"]) && is.finite(db["dom"]) &&
    da["dom"] != db["dom"] && da["pur"] >= 0.75 && db["pur"] >= 0.75

  # unknown-variance PATH via phi-remap (valid path test; the phi-remap UNION was
  # anti-conservative, removed).
  uv <- tryCatch(unknown_var_from_fit(fit, X, cfg$pair[1], cfg$pair[2]),
                 error = function(e) list(p_path_unknown = NA_real_))
  # VALID more-powerful unknown-variance UNION (and path) on the R-fiber.
  rf <- tryCatch(kmeans_union_unknownvar(X, k = cfg$k, cluster_1 = cfg$pair[1],
                   cluster_2 = cfg$pair[2], seed = cfg$init_seed, n_theta = 250, tol = 1e-7),
                 error = function(e) list(p_union = NA_real_, p_path = NA_real_))

  data.frame(
    rep            = r,
    p_union        = fit$pval_union,
    p_path         = fit$pval_path,
    p_naive        = fit$p_naive,
    p_path_unknown  = uv$p_path_unknown,
    p_union_rfib    = rf$p_union,
    p_path_rfib     = rf$p_path,
    test_stat      = fit$test_stat,
    n_paths        = length(fit$paths),
    exhaustive     = fit$exhaustive,
    detected       = unname(detected),
    runtime_s      = runtime,
    stringsAsFactors = FALSE)
}

# ---- cell setup (shared by raw + summary paths) ----------------------------
# Validates cfg and precomputes the per-cell constants (sizes, means, cov), so
# both run_cell() and run_cell_raw() -- and the checkpointed driver -- agree.
.cell_setup <- function(cfg) {
  cfg <- .fill_cfg(cfg)
  if (length(cfg$pair) != 2L || cfg$pair[1] == cfg$pair[2])
    stop("pair must be two distinct cluster labels")
  sizes <- .cluster_sizes(cfg)
  mu    <- .geometry_means(cfg$geometry, cfg$k, cfg$q, cfg$delta, cfg$elong)
  use_genCov <- !is.null(cfg$Sigma)
  SigChol <- NULL
  if (use_genCov) {
    if (!is.matrix(cfg$Sigma) || nrow(cfg$Sigma) != cfg$q || ncol(cfg$Sigma) != cfg$q)
      stop("Sigma must be a q x q matrix")
    SigChol <- chol(cfg$Sigma)             # X %*% SigChol has cov Sigma
    if (is.null(cfg$SigInv_test)) cfg$SigInv_test <- solve(cfg$Sigma)
  }
  list(cfg = cfg, sizes = sizes, n_tot = sum(sizes), mu = mu,
       use_genCov = use_genCov, SigChol = SigChol, is_null = .is_null_cell(cfg))
}

# Run a SPECIFIED set of replicate indices for a set-up cell; returns the raw
# per-rep data.frame (NULL-reps dropped). Exposed so the driver can run/resume
# reps in batches.  `rep_fun_wrap` lets the driver inject checkpointing.
.run_reps <- function(setup, idx, ncores = 1L, rep_fun_wrap = NULL) {
  f <- function(r) .one_rep(r, setup$cfg, setup$sizes, setup$mu,
                            setup$use_genCov, setup$SigChol)
  if (!is.null(rep_fun_wrap)) return(rep_fun_wrap(f, idx))   # driver-controlled
  reps <- mclapply(idx, f, mc.cores = ncores, mc.preschedule = FALSE)
  reps <- Filter(Negate(is.null), reps)
  if (length(reps) == 0) return(NULL)
  do.call(rbind, reps)
}

# Turn raw per-rep rows into the one-row cell summary.
.summarize_cell <- function(setup, res) {
  cfg <- setup$cfg; alpha <- cfg$alpha; is_null <- setup$is_null
  n_valid <- if (is.null(res)) 0L else nrow(res)
  rej <- function(p) if (n_valid == 0) NA_real_ else mean(p <= alpha, na.rm = TRUE)
  rej_union <- rej(res$p_union); rej_path <- rej(res$p_path); rej_naive <- rej(res$p_naive)
  rej_path_unk  <- rej(res$p_path_unknown)   # phi-remap path (valid)
  rej_union_rfib <- rej(res$p_union_rfib)    # R-fiber union (valid, more powerful)
  rej_path_rfib  <- rej(res$p_path_rfib)     # R-fiber path
  ks <- function(p) { p <- p[is.finite(p)]
    if (length(p) >= 2) suppressWarnings(ks.test(p, "punif")$p.value) else NA_real_ }
  ks_union <- NA_real_; ks_path <- NA_real_; ks_path_unk <- NA_real_
  ks_union_rfib <- NA_real_; ks_path_rfib <- NA_real_
  if (is_null && n_valid >= 2) {
    ks_union <- ks(res$p_union); ks_path <- ks(res$p_path)
    ks_path_unk <- ks(res$p_path_unknown)
    ks_union_rfib <- ks(res$p_union_rfib); ks_path_rfib <- ks(res$p_path_rfib)
  }

  # Chen & Witten estimands: detection probability (eq. 24) and CONDITIONAL
  # power (eq. 23) = rejection rate among reps where the tested pair was
  # recovered as two distinct true clusters.
  detected     <- if (!is.null(res$detected)) res$detected else logical(0)
  n_detected   <- sum(detected, na.rm = TRUE)
  detect_prob  <- if (n_valid == 0) NA_real_ else mean(detected, na.rm = TRUE)
  cond_pow <- function(p) {
    d <- detected & is.finite(p)
    if (sum(d) == 0) NA_real_ else mean(p[d] <= alpha)
  }
  cpow_union <- cond_pow(res$p_union); cpow_path <- cond_pow(res$p_path)
  cpow_path_unk <- cond_pow(res$p_path_unknown)
  cpow_union_rfib <- cond_pow(res$p_union_rfib); cpow_path_rfib <- cond_pow(res$p_path_rfib)
  data.frame(
    label            = if (is.null(cfg$label)) NA_character_ else cfg$label,
    geometry         = cfg$geometry,
    tested_pair_type = cfg$tested_pair_type,
    is_null          = is_null,
    k = cfg$k, q = cfg$q, delta = cfg$delta,
    n_total          = setup$n_tot,
    sizes            = paste(setup$sizes, collapse = "/"),
    pair             = paste(cfg$pair, collapse = "-"),
    cov_model        = if (setup$use_genCov) "genCov" else "iso",
    sig              = if (is.null(cfg$sig)) NA_real_ else cfg$sig,
    sig_test         = if (is.null(cfg$sig_test)) NA_real_ else cfg$sig_test,
    nreps            = cfg$nreps,
    n_valid          = n_valid,
    alpha            = alpha,
    rej_union        = rej_union, rej_path = rej_path, rej_naive = rej_naive,
    rej_path_unk     = rej_path_unk,
    rej_union_rfib   = rej_union_rfib, rej_path_rfib = rej_path_rfib,
    ks_union         = ks_union,  ks_path  = ks_path, ks_path_unk = ks_path_unk,
    ks_union_rfib    = ks_union_rfib, ks_path_rfib = ks_path_rfib,
    power_gain       = if (n_valid == 0) NA_real_ else (rej_union - rej_path),
    detect_prob      = detect_prob, n_detected = n_detected,
    cpow_union       = cpow_union, cpow_path = cpow_path, cpow_path_unk = cpow_path_unk,
    cpow_union_rfib  = cpow_union_rfib, cpow_path_rfib = cpow_path_rfib,
    cpow_gain        = if (is.na(cpow_union) || is.na(cpow_path)) NA_real_ else (cpow_union - cpow_path),
    mean_n_paths     = if (n_valid == 0) NA_real_ else mean(res$n_paths),
    frac_exhaustive  = if (n_valid == 0) NA_real_ else mean(res$exhaustive),
    median_runtime_s = if (n_valid == 0) NA_real_ else median(res$runtime_s),
    stringsAsFactors = FALSE)
}

# ---- run ONE cell: raw per-rep rows ---------------------------------------
# Returns the raw per-rep data.frame (with cfg fields not echoed; use the
# summary for those). Useful for pooling chunked sub-cells and for diagnostics.
run_cell_raw <- function(cfg) {
  setup <- .cell_setup(cfg)
  .run_reps(setup, seq_len(setup$cfg$nreps), ncores = setup$cfg$ncores)
}

# ---- run ONE cell: one-row summary (unchanged interface) -------------------
run_cell <- function(cfg) {
  setup <- .cell_setup(cfg)
  res   <- .run_reps(setup, seq_len(setup$cfg$nreps), ncores = setup$cfg$ncores)
  .summarize_cell(setup, res)
}

# ---- DRIVER: run a list of cells ------------------------------------------
# cells: a list of cfg lists. grid_ncores parallelizes ACROSS cells; each
# cell is run with its own (typically serial) ncores. Returns one rbind'd
# data.frame of per-cell summaries.
run_grid <- function(cells, grid_ncores = 1L, verbose = TRUE) {
  run_one <- function(i) {
    cfg <- cells[[i]]
    if (verbose) cat(sprintf("[cell %d/%d] %s\n", i, length(cells),
                             if (is.null(cfg$label)) "" else cfg$label))
    tryCatch(run_cell(cfg), error = function(e) {
      warning(sprintf("cell %d failed: %s", i, conditionMessage(e)))
      NULL
    })
  }
  out <- if (grid_ncores > 1L) {
    mclapply(seq_along(cells), run_one, mc.cores = grid_ncores,
             mc.preschedule = FALSE)
  } else {
    lapply(seq_along(cells), run_one)
  }
  out <- Filter(Negate(is.null), out)
  if (length(out) == 0) return(invisible(NULL))
  do.call(rbind, out)
}

# ---- script mode: tiny 2-cell smoke grid ----------------------------------
if (sys.nframe() == 0L && !interactive()) {
  cells <- list(
    list(label = "smoke-null",   geometry = "single-blob", delta = 0,
         n_per = 8L, k = 3L, q = 2L, pair = c(1, 3), sig = 1,
         tested_pair_type = "null-split", nreps = 30L,
         max_paths = 60L, ncores = 1L, cell_seed_offset = 0L),
    list(label = "smoke-strong", geometry = "triangle",   delta = 6,
         n_per = 8L, k = 3L, q = 2L, pair = c(1, 3), sig = 1,
         tested_pair_type = "different", nreps = 30L,
         max_paths = 60L, ncores = 1L, cell_seed_offset = 5000L)
  )
  t0 <- Sys.time()
  summ <- run_grid(cells, grid_ncores = 1L)
  cat(sprintf("\nDone in %.1fs\n",
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  cat("\n==== per-cell summary ====\n")
  print(summ, row.names = FALSE, digits = 3)
}
