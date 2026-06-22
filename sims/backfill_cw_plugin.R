#!/usr/bin/env Rscript
# ===========================================================================
# Backfill the Chen-Witten plug-in-variance p-values into the finished
# sweep_cw_uv chunks, following the SAME pattern used for p_union_rfib/p_path_rfib.
#
# Adds four columns to every (cell, rep) in sims/results/sweep_cw_uv/cells/*.rds:
#   p_union_med, p_path_med   (sigma_MED plug-in)
#   p_union_samp, p_path_samp (sigma_Sample plug-in)
# Does NOT modify any existing column. Idempotent. Atomic per-chunk writes.
#
# METHOD (important): the plug-in is the known-variance union/path test with
# sigma replaced by an estimate, on the SAME (sigma-independent) truncation set.
# We compute that set EXACTLY by re-running kmeans_inference_union (iso) and
# reusing fit$interval_union / fit$interval_path -- IDENTICAL to how the live
# harness .one_rep() now produces these columns (sweep_harness.R lines 263-272).
# This reproduces the stored p_union EXACTLY (max |diff| = 0 over 560 spot pairs),
# unlike the numerical phi-grid scan in kmeans_plugin_unknownvar(), whose 400-pt
# grid misses narrow intervals on ~13% of null reps (off by up to 0.7).
#
# Usage:
#   Rscript sims/backfill_cw_plugin.R verify [ncores]   # consistency check only
#   Rscript sims/backfill_cw_plugin.R run    [ncores]   # backfill all chunks
# ===========================================================================
suppressMessages(library(KmeansInference))
suppressMessages(library(parallel))
source("sims/sweep_harness.R")            # .cell_setup, .fill_cfg, .geometry_means, ...

CELLS_DIR <- "sims/results/sweep_cw_uv/cells"

# ---- reproduce build_grid("cw") WITHOUT running the sweep main block --------
load_build_grid <- function() {
  src <- parse("sims/sweep_run.R")
  env <- new.env()
  for (e in src) {
    if (is.call(e) && length(e) >= 1 && identical(e[[1]], as.name("<-"))) {
      lhs <- e[[2]]
      if (is.name(lhs) && as.character(lhs) %in% c("build_grid", ".chunks_of"))
        eval(e, env)
    }
  }
  env$build_grid
}

build_grid <- load_build_grid()
CELLS <- build_grid("cw")
LABELS_CFG <- setNames(CELLS, vapply(CELLS, function(c) c$label, character(1)))

# ---- regenerate X exactly as the harness does (.one_rep), returning the
#      FILLED cfg (sig_test = 1, etc.) ---------------------------------------
regen_X <- function(cfg, r) {
  setup <- .cell_setup(cfg); cfg <- setup$cfg
  n <- sum(setup$sizes)
  true_clusters <- rep(seq_len(cfg$k), times = setup$sizes)
  set.seed(cfg$data_seed0 + cfg$cell_seed_offset + r)
  Z <- matrix(rnorm(n * cfg$q), n, cfg$q)
  if (setup$use_genCov) {
    noise <- Z %*% setup$SigChol
  } else {
    sd_sim <- if (is.null(cfg$sig)) 1 else cfg$sig
    noise <- Z * sd_sim
  }
  X <- noise + setup$mu[true_clusters, , drop = FALSE]
  list(X = X, cfg = cfg)                 # cfg here is the FILLED cfg
}

# ---- EXACT plug-in p-values via the harness's own recipe -------------------
# Re-runs the iso union test, then evaluates the chi-ratio p-value on the EXACT
# interval_union / interval_path with sigma replaced by MED / Sample estimates.
# Returns p_union_known (== fit$pval_union) for the consistency check too.
plugin_exact <- function(X, fc) {
  fit <- suppressWarnings(suppressMessages(kmeans_inference_union(
    X, k = fc$k, cluster_1 = fc$pair[1], cluster_2 = fc$pair[2],
    iso = TRUE, sig = fc$sig_test, iter.max = fc$iter.max, seed = fc$init_seed,
    max_paths = fc$max_paths, scan_tail_prob = fc$scan_tail_prob, verbose = FALSE)))
  ivp <- function(iv, sf) tryCatch(
    KmeansInference:::.kmeans_interval_pvalue(iv, fit$test_stat, sf, fc$q),
    error = function(e) NA_real_)
  # ||nu||^2 = scale_factor / sig_test^2  (scale_factor = ||nu||^2 * sig_test^2)
  sqnu   <- fit$scale_factor / (if (is.null(fc$sig_test)) 1 else fc$sig_test)^2
  s_med  <- KmeansInference:::.kmeans_estimate_MED(X)
  s_samp <- sqrt(sum(scale(X, center = TRUE, scale = FALSE)^2) / (length(X) - fc$q))
  list(p_union_known = fit$pval_union,
       p_union_med  = ivp(fit$interval_union, sqnu * s_med^2),
       p_path_med   = ivp(fit$interval_path,  sqnu * s_med^2),
       p_union_samp = ivp(fit$interval_union, sqnu * s_samp^2),
       p_path_samp  = ivp(fit$interval_path,  sqnu * s_samp^2))
}

label_of_chunk <- function(fname) sub("__c[0-9]+\\.rds$", "", basename(fname))
chunk_files    <- function() list.files(CELLS_DIR, pattern = "__c[0-9]+\\.rds$",
                                        full.names = TRUE)

# ---- VERIFY: spot-check p_union_known vs stored p_union ---------------------
verify <- function(ncores = 1L) {
  pick <- c("cw_typeI_q2__c01.rds", "cw_typeI_q10__c01.rds",
            "cw_typeI_q50__c01.rds", "cw_typeI_q100__c01.rds",
            "cw_power_q10_d8__c01.rds", "cw_power_q50_d6__c01.rds",
            "cw_power_q2_d7__c01.rds")
  pick <- pick[file.exists(file.path(CELLS_DIR, pick))]
  jobs <- list()
  for (pf in pick) {
    ch <- readRDS(file.path(CELLS_DIR, pf))
    lab <- label_of_chunk(pf)
    for (r in head(ch$rep, 5L))
      jobs[[length(jobs) + 1L]] <- list(pf = pf, lab = lab, r = r,
                                        stored = ch$p_union[ch$rep == r])
  }
  f <- function(j) {
    cfg <- LABELS_CFG[[j$lab]]; gx <- regen_X(cfg, j$r); fc <- gx$cfg
    g <- tryCatch(plugin_exact(gx$X, fc), error = function(e) NULL)
    data.frame(chunk = j$pf, rep = j$r, q = fc$q,
               stored_p_union = j$stored,
               p_union_known = if (is.null(g)) NA_real_ else g$p_union_known,
               abs_diff = if (is.null(g)) NA_real_ else abs(j$stored - g$p_union_known),
               stringsAsFactors = FALSE)
  }
  res <- do.call(rbind, mclapply(jobs, f, mc.cores = ncores, mc.preschedule = FALSE))
  cat("=== consistency check: p_union_known vs stored p_union ===\n")
  print(res, row.names = FALSE, digits = 10)
  cat(sprintf("\nmax |abs_diff| = %.3e over %d pairs\n",
              max(res$abs_diff, na.rm = TRUE), nrow(res)))
  invisible(res)
}

# ---- one chunk: add the four columns (idempotent, atomic) ------------------
backfill_chunk <- function(fpath) {
  ch <- readRDS(fpath)
  if (nrow(ch) == 0) return(sprintf("%s (0 rows, skipped)", basename(fpath)))
  lab <- label_of_chunk(fpath)
  cfg <- LABELS_CFG[[lab]]
  if (is.null(cfg)) stop(sprintf("no cfg for label '%s'", lab))

  needed <- c("p_union_med", "p_path_med", "p_union_samp", "p_path_samp")
  if (all(needed %in% names(ch)) && !all(is.na(ch$p_union_med)))
    return(sprintf("%s (already done, skipped)", basename(fpath)))

  pum <- ppm <- pus <- pps <- rep(NA_real_, nrow(ch))
  for (i in seq_len(nrow(ch))) {
    gx <- regen_X(cfg, ch$rep[i]); fc <- gx$cfg
    g <- tryCatch(plugin_exact(gx$X, fc), error = function(e) NULL)
    if (!is.null(g)) {
      pum[i] <- g$p_union_med;  ppm[i] <- g$p_path_med
      pus[i] <- g$p_union_samp; pps[i] <- g$p_path_samp
    }
  }
  ch$p_union_med  <- pum; ch$p_path_med  <- ppm
  ch$p_union_samp <- pus; ch$p_path_samp <- pps

  tmp <- paste0(fpath, ".tmp-", Sys.getpid())
  saveRDS(ch, tmp); file.rename(tmp, fpath)
  sprintf("%s (%d rows, %d NA)", basename(fpath), nrow(ch), sum(is.na(pum)))
}

run_all <- function(ncores = 1L) {
  files <- chunk_files()
  cat(sprintf("[backfill] %d chunks, ncores=%d\n", length(files), ncores))
  t0 <- Sys.time()
  done <- mclapply(files, function(f) tryCatch(backfill_chunk(f),
            error = function(e) sprintf("%s ERROR: %s", basename(f), conditionMessage(e))),
          mc.cores = ncores, mc.preschedule = FALSE)
  done <- unlist(done)
  errs <- grep("ERROR", done, value = TRUE)
  cat(sprintf("[backfill] done %d chunks in %.1f min; %d errors\n",
              length(files), as.numeric(difftime(Sys.time(), t0, units = "mins")),
              length(errs)))
  if (length(errs)) { cat("ERRORS:\n"); cat(errs, sep = "\n"); cat("\n") }
  invisible(done)
}

# ---- main ------------------------------------------------------------------
if (sys.nframe() == 0L && !interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  mode <- if (length(args) >= 1) args[1] else "verify"
  ncores <- if (length(args) >= 2 && !is.na(as.integer(args[2]))) as.integer(args[2]) else 1L
  if (mode == "verify") verify(ncores)
  else if (mode == "run") run_all(ncores)
  else stop("mode must be 'verify' or 'run'")
}
