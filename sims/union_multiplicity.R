#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Empirical "effective multiplicity" of the union conditioning set.
#
# For each fit, the union set S_union = ⊔ I_j (disjoint phi-intervals). Under the
# null phi^2/tau ~ chi^2_q, so the reference mass of [a,b] is
#   pchisq(b^2/tau, q) - pchisq(a^2/tau, q).
# We report, per setting (averaged over reps):
#   n_paths : distinct Lloyd paths explored (raw, incl. negligible tail)
#   M_union : # disjoint intervals in S_union
#   M_eff   : participation ratio (sum w_j^2)^{-1}   <- "high-probability" count
#   M_ent   : exp(Shannon entropy of w)              <- alt effective count
#   M99     : min # regions covering 99% of union mass
#   top1    : mass fraction of the single largest region
#   eta     : share of EXTRA (union minus path) mass below phi_obs (power-useful)
#   waste   : 1 - M_eff/n_paths  (fraction of explored paths that are negligible)
# ---------------------------------------------------------------------------
suppressMessages({ library(KmeansInference); library(parallel) })

# reference chi mass of phi-interval matrix [a,b] given tau (scale_factor), df q.
# Cancellation-proof: take the max of the CDF-difference (stable in the left/center)
# and the SURVIVAL-difference (stable in the right tail, where both CDFs ~ 1).
imass <- function(mat, tau, q) {
  if (is.null(mat) || nrow(mat) == 0) return(numeric(0))
  lo_x <- mat[, 1]^2 / tau
  hi_x <- ifelse(is.infinite(mat[, 2]), Inf, mat[, 2]^2 / tau)
  m_cdf  <- pchisq(hi_x, q) - pchisq(lo_x, q)
  m_surv <- pchisq(lo_x, q, lower.tail = FALSE) - pchisq(hi_x, q, lower.tail = FALSE)
  pmax(m_cdf, m_surv, 0)
}

multiplicity <- function(fit) {
  tau <- fit$scale_factor; q <- length(fit$final_cluster) * 0 + NA  # set below
  q <- fit$df_q
  U  <- as.matrix(fit$interval_union)
  wj <- imass(U, tau, q); W <- sum(wj)
  if (W <= 0) return(NULL)
  w  <- wj / W
  M_eff <- 1 / sum(w^2)
  M_ent <- exp(-sum(ifelse(w > 0, w * log(w), 0)))
  sw <- sort(w, decreasing = TRUE)
  M99 <- which(cumsum(sw) >= 0.99)[1]
  # extra (union minus path) mass below phi_obs
  P  <- as.matrix(fit$interval_path)
  phi <- fit$test_stat
  below <- function(mat) { if (nrow(mat) == 0) return(0)
    m <- mat; m[, 2] <- pmin(m[, 2], phi); m <- m[m[, 1] < m[, 2], , drop = FALSE]
    sum(imass(m, tau, q)) }
  extra_total <- W - sum(imass(P, tau, q))
  eta <- if (extra_total > 1e-12) (below(U) - below(P)) / extra_total else NA_real_
  data.frame(n_paths = length(fit$paths), M_union = nrow(U),
             M_eff = M_eff, M_ent = M_ent, M99 = M99, top1 = max(w), eta = eta)
}

settings <- expand.grid(n_per = c(10L, 20L, 30L), q = c(2L, 10L), delta = c(2, 4))
nreps <- 30L; k <- 3L; sig <- 1
ncores <- if (length(commandArgs(TRUE)) >= 1) as.integer(commandArgs(TRUE)[1]) else 3L

one <- function(r, st) {
  set.seed(7000 + r)
  tc <- rep(1:k, each = st$n_per); n <- length(tc)
  mu <- rbind(c(st$delta/2,0), c(0, sqrt(3)*st$delta/2), c(-st$delta/2,0))
  M <- matrix(0, k, st$q); M[,1:2] <- mu
  X <- matrix(rnorm(n*st$q, sd = sig), n, st$q) + M[tc,]
  fit <- tryCatch(suppressWarnings(suppressMessages(
    kmeans_inference_union(X, k=k, cluster_1=1, cluster_2=3, sig=sig,
                           seed=2021, max_paths=200, verbose=FALSE))),
    error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  fit$df_q <- st$q
  multiplicity(fit)
}

out <- list()
for (i in seq_len(nrow(settings))) {
  st <- settings[i, ]
  rr <- mclapply(seq_len(nreps), one, st = st, mc.cores = ncores, mc.preschedule = FALSE)
  rr <- do.call(rbind, Filter(Negate(is.null), rr))
  if (is.null(rr) || nrow(rr) == 0) {
    cat(sprintf("setting n_per=%d q=%d delta=%g: 0 valid reps (skipped)\n",
                st$n_per, st$q, st$delta)); next
  }
  out[[i]] <- cbind(st, n_valid = nrow(rr),
                    as.data.frame(as.list(colMeans(as.matrix(rr), na.rm = TRUE))))
}
tab <- do.call(rbind, out)
tab$waste <- 1 - tab$M_eff / tab$n_paths
saveRDS(tab, "sims/results/union_multiplicity.rds")
cat("== effective multiplicity of the union (means over", nreps, "reps) ==\n")
print(tab[, c("n_per","q","delta","n_paths","M_union","M_eff","M99","top1","eta","waste")],
      row.names = FALSE, digits = 3)
