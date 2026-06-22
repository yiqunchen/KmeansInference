# Heavy-tailed robustness (Yun-Barber Figs 5-6): global-null Type I when the noise
# is t5 / t10 instead of Gaussian (variance standardized to 1). k=3, n_per=30, test
# clusters 1v2. Records known-sigma union/path + studentized R-fiber union/path + naive.
suppressMessages(library(KmeansInference)); source("sims/kmeans_union_unknownvar.R")
suppressMessages(library(parallel))
NREPS <- as.integer(Sys.getenv("NREPS","400")); NC <- 14
gen <- function(noise, n, q) {
  if (noise == "gaussian") matrix(rnorm(n*q), n, q)
  else { df <- as.integer(sub("t","",noise)); matrix(rt(n*q, df), n, q) * sqrt((df-2)/df) }
}
one <- function(rep, noise, q) {
  set.seed(7000 + 101*rep + 13*q)
  X <- gen(noise, 90, q)                      # mu = 0 : global null
  ft <- tryCatch(kmeans_inference_union(X, 3, 1, 2, sig = 1, seed = 2021), error = function(e) NULL)
  rf <- tryCatch(kmeans_union_unknownvar(X, 3, 1, 2, seed = 2021), error = function(e) NULL)
  if (is.null(ft) || is.null(rf)) return(NULL)
  data.frame(noise = noise, q = q, p_path = ft$pval_path, p_union = ft$pval_union,
             p_naive = ft$p_naive, p_path_rfib = rf$p_path, p_union_rfib = rf$p_union)
}
rows <- list()
for (noise in c("gaussian","t5","t10")) for (q in c(2,10)) {
  r <- mclapply(1:NREPS, function(i) tryCatch(one(i, noise, q), error=function(e) NULL), mc.cores = NC)
  r <- do.call(rbind, Filter(Negate(is.null), r)); rows[[paste(noise,q)]] <- r
  ti <- function(p) mean(p <= 0.05, na.rm=TRUE)
  cat(sprintf("%-9s q=%2d n=%3d | known: path=%.3f union=%.3f | studF: path=%.3f union=%.3f | naive=%.3f\n",
      noise, q, nrow(r), ti(r$p_path), ti(r$p_union), ti(r$p_path_rfib), ti(r$p_union_rfib), ti(r$p_naive)))
  flush.console()
}
saveRDS(do.call(rbind, rows), "sims/results/heavytail_typeI.rds")
cat("DONE -> sims/results/heavytail_typeI.rds\n")
