# ---------------------------------------------------------------------------
# plot_exact_vs_is.R -- our EXACT union p-value vs Yun-Barber IMPORTANCE SAMPLING.
#
# Yun & Barber (2023) importance-sample the K>2 truncation set (their
# fun_proposed_approx, ndraws=8000 in their paper) because they lack an explicit
# characterization. Our R-fiber solver computes the same truncated-F probability
# EXACTLY. Here we run THEIR importance sampler verbatim -- proposal, R-fiber
# rotation, importance weights -- swapping only the clustering oracle from HAC to
# k-means (and to our two-tested-cluster preservation), so it estimates exactly
# our p_union. Repeated at several draw counts, the IS estimate is noisy and
# converges to our exact value; our value is the zero-variance limit.
# ---------------------------------------------------------------------------
source("sims/house_style.R")
suppressMessages({ library(ggplot2); library(truncnorm); library(intervals); library(parallel) })
source("/tmp/yjyun_repo/functions.R")          # fun_P0/P1/P2, fun_ts, fun_ratio_target_to_proposal (Yun-Barber)
suppressMessages(library(KmeansInference))
source("sims/unknown_var_union.R"); source("sims/rfiber_exact.R")
source("sims/kmeans_union_unknownvar_exact.R") # .build_fiber, .lloyd_, .presv_, exact solver

SEED <- 17; D <- 3.0; NPER <- 14; Q <- 2; K <- 3; C1 <- 1; C2 <- 2; ALPHA <- 0.05
set.seed(SEED)
mu <- rbind(c(0,0), c(D,0), c(D/2, D*sqrt(3)/2))
X  <- do.call(rbind, lapply(1:3, function(g) matrix(rnorm(NPER*Q),NPER,Q) + matrix(mu[g,],NPER,Q,byrow=TRUE)))
ex <- kmeans_union_unknownvar_exact(X, K, C1, C2, seed = SEED)
p_exact <- ex$p_union
est <- kmeans_estimation(X, K, 10, SEED, verbose = FALSE); cl <- est$final_cluster; init <- est$random_init_obs

# Yun-Barber's importance sampler, clustering oracle swapped to k-means + our
# two-tested-cluster preservation (matches the exact solver's union set).
is_kmeans <- function(ndraws) {
  q <- ncol(X); n1 <- sum(cl==C1); n2 <- sum(cl==C2); m <- n1+n2
  P0 <- fun_P0(cl,C1,C2); P1 <- fun_P1(cl,C1,C2); P2 <- fun_P2(cl,C1,C2)
  psi0 <- sqrt(sum((P0%*%X)^2)); psi1 <- sqrt(sum((P1%*%X)^2))
  ts <- fun_ts(X, cl, C1, C2); ts_Beta <- (ts/(m-2))/(1+ts/(m-2))
  samp <- rtruncnorm(ndraws, a=0, b=1, mean=ts_Beta, sd=ALPHA)
  w <- sapply(samp, fun_ratio_target_to_proposal, ts_=ts_Beta, m=m, q=q, alpha=ALPHA)
  vec_h1 <- as.integer(samp >= ts_Beta)
  t1 <- sqrt(psi0^2+psi1^2); t2 <- P0%*%X/psi0; t3 <- P1%*%X/psi1
  vec_h2 <- integer(ndraws)
  for (j in 1:ndraws) {
    y <- t1*sqrt(samp[j])*t2 + t1*sqrt(1-samp[j])*t3 + P2%*%X
    cy <- tryCatch(.lloyd_(y, K, init, 10, tol_eps=1e-6, verbose=FALSE)$final_cluster, error=function(e) NULL)
    if (!is.null(cy) && .presv_(cl, cy, C1, C2)) vec_h2[j] <- 1L
  }
  vec_h1[vec_h2==0] <- 0L
  w_new <- w[vec_h2==1]; w_new <- w_new - max(w_new)
  sum(exp(w_new[vec_h1[vec_h2==1]==1]) / sum(exp(w_new)))
}

NDRAWS <- c(250, 500, 1000, 2000, 4000, 8000); REPS <- 20
jobs <- expand.grid(nd = NDRAWS, rep = 1:REPS)
nc <- max(1L, min(12L, detectCores() - 2L))
cat(sprintf("exact p_union=%.4f ; running IS (%d draw-counts x %d reps) on %d cores\n", p_exact, length(NDRAWS), REPS, nc))
ests <- unlist(mclapply(seq_len(nrow(jobs)), function(i) tryCatch(is_kmeans(jobs$nd[i]), error=function(e) NA_real_), mc.cores = nc))
df <- data.frame(ndraws = factor(jobs$nd), p = ests)
saveRDS(list(df = df, p_exact = p_exact), "sims/results/exact_vs_is.rds")

fig <- ggplot(df, aes(ndraws, p)) +
  geom_hline(yintercept = p_exact, colour = "black", linewidth = 0.9) +
  geom_boxplot(fill = km_col[["studentized"]], alpha = 0.25, colour = km_col[["studentized"]],
               outlier.size = 0.5, width = 0.6, linewidth = 0.45) +
  annotate("text", x = 0.7, y = p_exact, vjust = -0.7, hjust = 0, fontface = "bold", size = 3.4,
           label = sprintf("our exact  p = %.4f", p_exact)) +
  labs(x = "importance-sampling draws (Yun--Barber)", y = "estimated union p-value",
       title = "Exact union p-value vs. importance sampling (K=3)") +
  theme_km()
ggsave_km(fig, "sims/results/exact_vs_is", width = 7.4, height = 4.6)
cat("Wrote exact_vs_is.{pdf,png}\n")
