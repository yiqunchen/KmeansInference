# ---------------------------------------------------------------------------
# validate_exact_backfill.R -- compare the EXACT R-fiber re-run against the
# NUMERICAL run, cell by cell.
#
# Determinism check: the known-variance + plug-in columns are solver-independent,
# so with identical seeds they MUST match the numerical run to ~machine precision.
# If they do, the only thing the backfill changed is the R-fiber (*_rfib) columns,
# and it is safe to promote sweep_cw_uv_exact -> canonical.
#
# Usage: Rscript sims/validate_exact_backfill.R
# ---------------------------------------------------------------------------
NUM <- "sims/results/sweep_cw_uv/summary.rds"
EXA <- "sims/results/sweep_cw_uv_exact/summary.rds"
stopifnot(file.exists(NUM), file.exists(EXA))

num <- as.data.frame(readRDS(NUM))
exa <- as.data.frame(readRDS(EXA))
key <- "label"
common <- intersect(num[[key]], exa[[key]])
cat(sprintf("cells: numerical %d, exact %d, common %d\n", nrow(num), nrow(exa), length(common)))
num <- num[match(common, num[[key]]), ]; exa <- exa[match(common, exa[[key]]), ]

# columns that MUST be identical (solver-independent)
known_cols <- intersect(names(num), c(
  "rej_naive", "rej_path", "rej_union", "rej_path_unk",
  "ks_path", "ks_union", "cpow_path", "cpow_union",
  "rej_union_med", "rej_path_med", "rej_union_samp", "rej_path_samp",
  "detect_prob", "power_gain"))
# columns EXPECTED to move (the backfill target)
rfib_cols <- intersect(names(num), c(
  "rej_union_rfib", "rej_path_rfib", "ks_union_rfib", "ks_path_rfib",
  "cpow_union_rfib", "cpow_path_rfib"))

maxdiff <- function(cols) sapply(cols, function(c) {
  a <- suppressWarnings(as.numeric(num[[c]])); b <- suppressWarnings(as.numeric(exa[[c]]))
  ok <- is.finite(a) & is.finite(b)
  if (!any(ok)) return(NA_real_)
  max(abs(a[ok] - b[ok]))
})

cat("\n--- KNOWN / plug-in columns (must match; expect ~0) ---\n")
kd <- maxdiff(known_cols); print(round(kd, 8))
cat("\n--- R-fiber columns (backfill target; expect to move) ---\n")
rd <- maxdiff(rfib_cols); print(round(rd, 6))

tol <- 1e-6
bad <- known_cols[is.finite(kd) & kd > tol]
cat("\n=== VERDICT ===\n")
if (length(bad) == 0) {
  cat(sprintf("PASS: all known/plug-in columns match within %g (max %g).\n",
              tol, max(kd, na.rm = TRUE)))
  cat("The exact run reproduces the numerical run except on the R-fiber columns;\n")
  cat("safe to promote sweep_cw_uv_exact -> canonical.\n")
} else {
  cat("FAIL: these known columns differ beyond tolerance -- investigate before promoting:\n")
  for (c in bad) cat(sprintf("  %s: max|diff| = %g\n", c, kd[[c]]))
}
