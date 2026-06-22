# Reproduce Yun-Barber (2301.12999) fig_power_K3_equidistant SETTINGS in our
# k-means framework: n_each=10, K=3, equilateral means, sigma=1, delta 0..7.
# Q (dimension) via env var Q (default 2); reps via NREPS (default 300).
#
# APPLES-TO-APPLES variance comparison (Yun-Barber's design): hold the truncation
# set FIXED, vary ONLY the sigma-handling.  All PATH methods use the SAME phi-line
# path set fit$interval_path:
#   or_p = oracle (true sigma, chi pivot)
#   uf_p = studentized-F on the SAME set (rej_path_unk; unknown sigma, F pivot)  <- fair "unknown" curve
#   md_p = sigma_MED plug-in (chi pivot)
#   sa_p = sigma_all (sample) plug-in (chi pivot)
# On a fixed set, oracle >= studentized-F >= plug-ins (known sigma is the best).
#
# CROSS-CONDITIONING (to explain why R-fiber F looked like it beat the oracle):
#   rf_p = R-fiber path F (DIFFERENT, weaker conditioning -- can exceed phi oracle)
# UNION (conditioning axis): or_u/md_u/sa_u are phi-line union (known sigma);
#   rf_u = R-fiber union F (the only VALID unknown-sigma union).
source("sims/sweep_harness.R"); suppressMessages(library(parallel))
DELTAS <- 0:7; NREPS <- as.integer(Sys.getenv("NREPS", "300")); Q <- as.integer(Sys.getenv("Q", "2"))
base <- list(k = 3, q = Q, pair = c(1, 2), sig = 1, n_per = 10, max_paths = 250, ncores = 7)
rows <- list()
for (d in DELTAS) {
  r <- run_cell(c(base, list(geometry = "triangle", delta = d, nreps = NREPS)))
  rows[[length(rows)+1]] <- data.frame(delta = d, n_valid = r$n_valid, detect = r$detect_prob,
    # PATH (same phi-line set; vary sigma-handling) + R-fiber cross-conditioning
    or_p = r$rej_path, uf_p = r$rej_path_unk, md_p = r$rej_path_med, sa_p = r$rej_path_samp,
    rf_p = r$rej_path_rfib,
    # UNION (phi-line known-sigma) + valid R-fiber F-union
    or_u = r$rej_union, md_u = r$rej_union_med, sa_u = r$rej_union_samp, rf_u = r$rej_union_rfib)
  cat(sprintf("q=%d d=%d n=%3d det=%.2f | PATH(sameS) oracle=%.3f studF=%.3f MED=%.3f samp=%.3f [Rfib=%.3f] | UNION oracle=%.3f Rfib=%.3f\n",
      Q, d, r$n_valid, r$detect_prob, r$rej_path, r$rej_path_unk, r$rej_path_med, r$rej_path_samp,
      r$rej_path_rfib, r$rej_union, r$rej_union_rfib)); flush.console()
}
res <- do.call(rbind, rows)
saveRDS(res, sprintf("sims/results/repro_yunbarber_q%d.rds", Q))
cat(sprintf("\nSaved sims/results/repro_yunbarber_q%d.rds\n", Q))
