# ---------------------------------------------------------------------------
# checkpoint.R  --  robust, interrupt-safe incremental saving for simulations.
#
# Design goals
#   * Every completed unit of work is flushed to disk immediately, so a killed
#     run (Ctrl-C, crash, OOM, wall-clock limit) never loses finished work.
#   * Writes are ATOMIC (write tempfile in the same dir, then file.rename), so a
#     process dying mid-write can never leave a half-written / corrupt file.
#   * Resume is automatic: re-running skips any unit whose checkpoint exists.
#   * No external dependencies (base R only); one file per unit so parallel
#     workers (mclapply forks) never clobber each other.
#
# Two granularities are supported:
#   - CELL level  : one summary row per parameter cell  -> cp$cells_dir/<key>.rds
#   - REP  level  : raw per-replicate rows within a cell -> cp$raw_dir/<key>.rds
#                   (so a single long cell also resumes mid-way, batch by batch)
# ---------------------------------------------------------------------------

# Open (or re-open) a checkpoint store rooted at `dir`.
cp_open <- function(dir) {
  cells_dir <- file.path(dir, "cells")
  raw_dir   <- file.path(dir, "raw")
  dir.create(cells_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(raw_dir,   recursive = TRUE, showWarnings = FALSE)
  structure(list(dir = dir, cells_dir = cells_dir, raw_dir = raw_dir),
            class = "cp_store")
}

# A filesystem-safe key. Use a short, stable label for each cell.
cp_key <- function(label) gsub("[^A-Za-z0-9._-]", "_", label)

# Atomic single-object save: temp file in the same directory, then rename.
.cp_atomic_saveRDS <- function(obj, path) {
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(obj, tmp)
  file.rename(tmp, path)   # atomic on the same filesystem
  invisible(path)
}

# ---- CELL level -------------------------------------------------------------
cp_cell_path <- function(cp, key) file.path(cp$cells_dir, paste0(cp_key(key), ".rds"))
cp_has_cell  <- function(cp, key) file.exists(cp_cell_path(cp, key))
cp_save_cell <- function(cp, key, row) .cp_atomic_saveRDS(row, cp_cell_path(cp, key))

# Collect every saved cell summary into one data.frame (NULL if none yet).
cp_collect_cells <- function(cp) {
  files <- list.files(cp$cells_dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(files) == 0) return(NULL)
  rows <- lapply(files, function(f) tryCatch(readRDS(f), error = function(e) NULL))
  # keep only non-empty frames; align on the full (modal-max) column schema so a
  # 0-valid-rep chunk saved as a stub frame can never break the rbind.
  rows <- Filter(function(x) !is.null(x) && nrow(x) > 0, rows)
  if (length(rows) == 0) return(NULL)
  nc <- vapply(rows, ncol, 1L)
  rows <- rows[nc == max(nc)]
  do.call(rbind, rows)
}

# ---- REP level (resume a single cell batch-by-batch) ------------------------
cp_raw_path <- function(cp, key) file.path(cp$raw_dir, paste0(cp_key(key), ".rds"))

# Run `nreps` replicates of `rep_fun(i)` for one cell, flushing accumulated raw
# rows to disk after every `batch` reps. On resume, already-saved reps are read
# back and skipped, so only the remaining reps run. `rep_fun(i)` must return a
# 1-row data.frame (or NULL to drop that rep). Returns all raw rows rbind-ed.
cp_run_reps <- function(cp, key, nreps, rep_fun, batch = 25L, log_fun = NULL) {
  path <- cp_raw_path(cp, key)
  done <- if (file.exists(path)) tryCatch(readRDS(path), error = function(e) NULL) else NULL
  done_ids <- if (!is.null(done) && "rep" %in% names(done)) unique(done$rep) else integer(0)
  todo <- setdiff(seq_len(nreps), done_ids)
  acc  <- if (is.null(done)) list() else list(done)
  if (length(todo) == 0) return(do.call(rbind, acc))

  chunks <- split(todo, ceiling(seq_along(todo) / batch))
  for (ch in chunks) {
    new_rows <- lapply(ch, function(i) tryCatch(rep_fun(i), error = function(e) NULL))
    new_rows <- Filter(Negate(is.null), new_rows)
    if (length(new_rows) > 0) acc <- c(acc, list(do.call(rbind, new_rows)))
    combined <- do.call(rbind, acc)
    .cp_atomic_saveRDS(combined, path)        # flush after each batch
    if (!is.null(log_fun))
      log_fun(sprintf("cell %s: %d/%d reps", key, nrow(combined), nreps))
  }
  do.call(rbind, acc)
}

# Rebuild and persist the combined results table from all saved cells.
# Writes <dir>/summary.rds and <dir>/summary.csv atomically; returns the table.
cp_finalize <- function(cp, csv = TRUE) {
  tab <- cp_collect_cells(cp)
  if (is.null(tab)) return(invisible(NULL))
  .cp_atomic_saveRDS(tab, file.path(cp$dir, "summary.rds"))
  if (csv) {
    tmp <- file.path(cp$dir, paste0("summary.csv.tmp-", Sys.getpid()))
    utils::write.csv(tab, tmp, row.names = FALSE)
    file.rename(tmp, file.path(cp$dir, "summary.csv"))
  }
  tab
}
