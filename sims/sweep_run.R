#!/usr/bin/env Rscript
# ===========================================================================
# Comprehensive, INTERRUPT-SAFE sweep driver: union vs path vs naive.
#
# Mirrors Chen & Witten (2022, arXiv:2203.15267) Section 5 design and adds our
# extensions (general covariance, unbalanced sizes, k-sweep). The union test's
# pval_union is compared against pval_path (= their p_selective) and p_naive.
#
# ROBUSTNESS / SAFE INTERRUPT (sims/checkpoint.R):
#   * Each logical cell's nreps is split into CHUNKS of `chunk` reps. Every
#     chunk is a checkpoint unit: its raw per-rep rows are flushed atomically to
#     <ckpt>/cells/<label>__cNN.rds the moment it finishes.
#   * A killed run loses at most one in-flight chunk per core; re-running SKIPS
#     every chunk already on disk and continues. Same script + same grid => exact
#     resume. Chunks run in parallel across cores; each writes its own file.
#   * finalize() pools the raw rows per logical cell and writes
#     <ckpt>/summary.{rds,csv} (atomic). Safe to run any time, even mid-sweep.
#
# Usage:
#   Rscript sims/sweep_run.R [tier] [grid_ncores] [chunk] [ckpt_dir]
#     tier        : smoke | standard | full          (default standard)
#     grid_ncores : parallel chunks                  (default cores-1)
#     chunk       : reps per checkpoint unit         (default 50)
#     ckpt_dir    : checkpoint root                  (default sims/results/sweep_<tier>)
#   # resume: just re-run the SAME command. finalize-only: append 'finalize'.
# ===========================================================================
source("sims/checkpoint.R")
source("sims/sweep_harness.R")

# ---- the comprehensive grid (list of logical cells) ------------------------
# Each cell is a cfg list with a UNIQUE `label`. nreps & a few sizes scale with
# tier. Geometry/q/delta/sigma/k follow Chen & Witten Sec. 5 where noted.
build_grid <- function(tier) {
  ## --- Chen & Witten (2022) Section 5, expanded: n=90 (30/cluster), K=3.
  ## Type I (Fig 4): global null mu=0, q in {2,10,50,100}, 1000 reps (tight QQ/CI).
  ## Power  (Fig 5): triangle, q in {2,10,50}, sigma=1 (non-saturating),
  ##   integer delta=1..8 (fine where power climbs), 200 reps (SE bars in plot).
  ## Reports detection probability + conditional power (eqs. 23-24) and known +
  ## unknown variance for every method.
  if (tier == "cw") {
    cells <- list(); off <- 0L
    addcw <- function(cfg) { off <<- off + 1000L
      cfg$cell_seed_offset <- off; cfg$n_per <- 30L; cfg$k <- 3L
      cfg$max_paths <- 400L; cfg$ncores <- 1L; cfg$pair <- c(1, 3)
      cells[[length(cells) + 1L]] <<- cfg }
    for (q in c(2L, 10L, 50L, 100L))
      addcw(list(label = sprintf("cw_typeI_q%d", q), geometry = "single-blob",
                 delta = 0, q = q, sig = 1, tested_pair_type = "null-split",
                 nreps = 1000L))
    for (q in c(2L, 10L, 50L)) for (dl in 1:8)
      addcw(list(label = sprintf("cw_power_q%d_d%d", q, dl), geometry = "triangle",
                 delta = dl, q = q, sig = 1, tested_pair_type = "different",
                 nreps = 200L))
    return(cells)
  }

  reps_typeI <- switch(tier, smoke = 40L, standard = 400L, full = 2000L)
  reps_power <- switch(tier, smoke = 40L, standard = 200L, full = 1000L)
  n_per      <- switch(tier, smoke = 8L,  standard = 20L,  full = 30L)
  mp         <- switch(tier, smoke = 60L, standard = 150L, full = 200L)
  cells <- list(); off <- 0L
  add <- function(cfg) { off <<- off + 1000L
    cfg$cell_seed_offset <- off; cfg$max_paths <- mp; cfg$ncores <- 1L
    cells[[length(cells) + 1L]] <<- cfg }

  ## --- Panel A: selective Type I error (global null), vs dimension q -------
  ## C&W Fig. 4: mu = 0, K=3, sig=1, q in {2,10,50,100}. (q capped by n_per*k.)
  for (q in c(2L, 10L, 50L)) {
    add(list(label = sprintf("typeI_q%d", q), geometry = "single-blob",
             delta = 0, k = 3L, q = q, n_per = n_per, sig = 1,
             pair = c(1, 3), tested_pair_type = "null-split", nreps = reps_typeI))
  }

  ## --- Panel B: power vs separation delta, by sigma (C&W Fig. 5) -----------
  ## triangle, K=3, q=10 (and q=2), delta=2..8, sigma in {0.5,1}.
  for (q in c(2L, 10L)) for (sg in c(0.5, 1)) for (dl in c(2, 3, 4, 5, 6, 8)) {
    add(list(label = sprintf("power_q%d_s%g_d%g", q, sg, dl),
             geometry = "triangle", delta = dl, k = 3L, q = q, n_per = n_per,
             sig = sg, pair = c(1, 3), tested_pair_type = "different",
             nreps = reps_power))
  }

  ## --- Panel C: number of clusters k (our extension) ----------------------
  for (kk in c(2L, 4L, 5L)) {
    pr <- if (kk == 2L) c(1, 2) else c(1, 3)
    add(list(label = sprintf("k%d", kk), geometry = "triangle", delta = 5,
             k = kk, q = 2L, n_per = n_per, sig = 1, pair = pr,
             tested_pair_type = "different", nreps = reps_power))
  }

  ## --- Panel D: general covariance (our extension) ------------------------
  ## AR(1) and heteroscedastic Sigma; include delta=0 for genCov Type I.
  ar1 <- function(rho, q) rho^abs(outer(seq_len(q), seq_len(q), "-"))
  for (rho in c(0.5, 0.9)) for (dl in c(0, 5)) {
    add(list(label = sprintf("gencov_ar%g_d%g", rho, dl), geometry = "triangle",
             delta = dl, k = 3L, q = 2L, n_per = n_per, Sigma = ar1(rho, 2L),
             pair = c(1, 3),
             tested_pair_type = if (dl == 0) "null-split" else "different",
             nreps = if (dl == 0) reps_typeI else reps_power))
  }

  ## --- Panel E: unbalanced cluster sizes (our extension) ------------------
  for (dl in c(0, 5)) {
    add(list(label = sprintf("unbal_d%g", dl), geometry = "triangle", delta = dl,
             k = 3L, q = 2L, n = 3L * n_per, balance = c(1, 2, 3), sig = 1,
             pair = c(1, 3),
             tested_pair_type = if (dl == 0) "null-split" else "different",
             nreps = if (dl == 0) reps_typeI else reps_power))
  }
  cells
}

# ---- expand a logical cell into checkpointed chunks ------------------------
.chunks_of <- function(nreps, chunk) {
  if (nreps <= chunk) return(list(seq_len(nreps)))
  split(seq_len(nreps), ceiling(seq_len(nreps) / chunk))
}

# ---- main ------------------------------------------------------------------
args        <- commandArgs(trailingOnly = TRUE)
tier        <- if (length(args) >= 1) args[1] else "standard"
finalize_only <- "finalize" %in% args
grid_ncores <- if (length(args) >= 2 && !is.na(as.integer(args[2]))) as.integer(args[2]) else max(1L, detectCores() - 1L)
chunk       <- if (length(args) >= 3 && !is.na(as.integer(args[3]))) as.integer(args[3]) else 50L
ckpt_dir    <- if (length(args) >= 4) args[4] else file.path("sims/results", paste0("sweep_", tier))

cells <- build_grid(tier)
cp    <- cp_open(ckpt_dir)
labels_cfg <- setNames(cells, vapply(cells, function(c) c$label, character(1)))

# Build the flat work list of (label, chunk-index, rep-indices), skipping done.
work <- list()
for (cfg in cells) {
  setup <- .cell_setup(cfg)           # validate early so a bad cell fails fast
  chs <- .chunks_of(setup$cfg$nreps, chunk)
  for (ci in seq_along(chs)) {
    key <- sprintf("%s__c%02d", cfg$label, ci)
    if (!cp_has_cell(cp, key))
      work[[length(work) + 1L]] <- list(key = key, label = cfg$label,
                                        cfg = cfg, idx = chs[[ci]])
  }
}

cat(sprintf("[sweep_run] tier=%s | %d logical cells | %d chunks to run (chunk=%d) | ncores=%d | ckpt=%s\n",
            tier, length(cells), length(work), chunk, grid_ncores, ckpt_dir))

if (!finalize_only && length(work) > 0) {
  run_chunk <- function(w) {
    setup <- .cell_setup(w$cfg)
    raw <- .run_reps(setup, w$idx, ncores = 1L)     # serial within a chunk
    if (is.null(raw) || nrow(raw) == 0) {
      raw <- data.frame(label = character(0))       # clean 0-row marker
    } else {
      raw$label <- w$label
    }
    cp_save_cell(cp, w$key, raw)                     # atomic flush
    sprintf("%s (%d reps)", w$key, if (is.null(raw)) 0L else nrow(raw))
  }
  done <- if (grid_ncores > 1L)
    mclapply(work, run_chunk, mc.cores = grid_ncores, mc.preschedule = FALSE)
  else lapply(work, run_chunk)
  cat(sprintf("[sweep_run] ran %d chunks\n", length(done)))
}

# ---- finalize: pool raw rows per logical cell -> summary -------------------
all_raw <- cp_collect_cells(cp)
if (is.null(all_raw) || nrow(all_raw) == 0) { cat("[sweep_run] no results yet\n"); quit(save = "no") }

summ_rows <- list()
for (lab in unique(all_raw$label)) {
  setup <- .cell_setup(labels_cfg[[lab]])
  res   <- all_raw[all_raw$label == lab, , drop = FALSE]
  res   <- res[, setdiff(names(res), "label"), drop = FALSE]
  if (nrow(res) == 0) next
  summ_rows[[lab]] <- .summarize_cell(setup, res)
}
summary_tab <- do.call(rbind, summ_rows)
.cp_atomic_saveRDS(summary_tab, file.path(ckpt_dir, "summary.rds"))
tmp <- file.path(ckpt_dir, paste0("summary.csv.tmp-", Sys.getpid()))
utils::write.csv(summary_tab, tmp, row.names = FALSE)
file.rename(tmp, file.path(ckpt_dir, "summary.csv"))

cat(sprintf("\n[sweep_run] finalized %d cells -> %s/summary.csv\n",
            nrow(summary_tab), ckpt_dir))
print(summary_tab[, c("label", "is_null", "k", "q", "delta", "n_valid",
                      "rej_union", "rej_path", "rej_naive", "power_gain",
                      "ks_union", "mean_n_paths")],
      row.names = FALSE, digits = 3)
