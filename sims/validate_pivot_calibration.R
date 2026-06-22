# DIRECT end-to-end R-fiber pivot uniformity at the GLOBAL NULL (mu=0 => H0' holds
# exactly => p MUST be ~Uniform if the solver is correct). q=2 vs q=10, no union
# exploration (isolates the R-fiber). KS + Type I at several alpha.
suppressMessages({library(KmeansInference); library(parallel)})
source("sims/kmeans_union_unknownvar.R")
NREPS <- as.integer(Sys.getenv("NREPS","1500"))
for (q in c(2, 10)) {
  res <- mclapply(1:NREPS, function(i) {
    set.seed(60000 + i + 7L*q)
    X <- matrix(rnorm(90*q), 90, q)                 # global null
    rf <- tryCatch(kmeans_union_unknownvar(X, 3, 1, 2, seed = 2021, n_theta = 300), error=function(e) NULL)
    if (is.null(rf)) return(NULL)
    c(path = rf$p_path, union = rf$p_union, R_obs = rf$R_obs, nint = rf$n_intervals)
  }, mc.cores = 14)
  m <- do.call(rbind, Filter(Negate(is.null), res))
  for (nm in c("path","union")) { p <- m[,nm]; p <- p[is.finite(p)]
    cat(sprintf("q=%2d Rfib_%-5s | a.01=%.3f a.05=%.3f a.10=%.3f a.20=%.3f a.50=%.3f | KSp=%.3f  (n=%d, med_nint=%.1f, med_Robs=%.1f)\n",
      q, nm, mean(p<=.01), mean(p<=.05), mean(p<=.10), mean(p<=.20), mean(p<=.50),
      suppressWarnings(ks.test(p[p>0&p<1]+1e-9*rnorm(sum(p>0&p<1)),"punif")$p.value), length(p), median(m[,"nint"],na.rm=TRUE), median(m[,"R_obs"],na.rm=TRUE)))
    flush.console() }
}
cat("DONE\n")
