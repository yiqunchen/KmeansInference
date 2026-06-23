# ---------------------------------------------------------------------------
# power_geometries.R -- known-sigma union vs path power across cluster GEOMETRIES,
# following Yun & Barber (2023) Fig 4: K=2 horizontal, K=3 equilateral triangle,
# K=3 collinear. Shows the union gain is robust to cluster configuration.
# Saves sims/results/power_geometries.rds ; plot with plot_power_geometries.R.
# ---------------------------------------------------------------------------
suppressMessages({ library(KmeansInference); library(parallel) })

NPER <- 20L; Q <- 2L; REPS <- 150L; DELTAS <- 0:7
geoms <- list(
  "K=2 horizontal" = function(d) list(mu = rbind(c(0,0), c(d,0)),              k = 2L, pair = c(1,2)),
  "K=3 triangle"   = function(d) list(mu = rbind(c(0,0), c(d,0), c(d/2, d*sqrt(3)/2)), k = 3L, pair = c(1,2)),
  "K=3 collinear"  = function(d) list(mu = rbind(c(0,0), c(d,0), c(2*d,0)),    k = 3L, pair = c(1,2)))

gen <- function(sp, seed) { set.seed(seed); K <- nrow(sp$mu)
  do.call(rbind, lapply(1:K, function(j) matrix(rnorm(NPER*Q), NPER, Q) + matrix(sp$mu[j,], NPER, Q, byrow = TRUE))) }

grid <- expand.grid(g = names(geoms), d = DELTAS, rep = 1:REPS, stringsAsFactors = FALSE)
gi   <- setNames(seq_along(geoms), names(geoms))
run1 <- function(i) {
  gname <- grid$g[i]; d <- grid$d[i]; rp <- grid$rep[i]
  sp <- geoms[[gname]](d); seed <- gi[[gname]]*1e6 + d*1e3 + rp
  X  <- gen(sp, seed)
  fit <- tryCatch(kmeans_inference_union(X, sp$k, sp$pair[1], sp$pair[2], sig = 1, seed = seed),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  data.frame(geom = gname, delta = d, rep = rp, p_path = fit$pval_path, p_union = fit$pval_union)
}
nc <- max(1L, min(12L, detectCores() - 2L))
cat(sprintf("[geom] %d cells x %d reps over %d geometries on %d cores\n", length(DELTAS), REPS, length(geoms), nc))
res <- mclapply(seq_len(nrow(grid)), run1, mc.cores = nc)
df  <- do.call(rbind, res[!vapply(res, is.null, logical(1))])
saveRDS(df, "sims/results/power_geometries.rds")
cat(sprintf("[geom] done: %d rows -> sims/results/power_geometries.rds\n", nrow(df)))
