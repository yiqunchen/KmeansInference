#!/usr/bin/env Rscript
# Test the assembler (Module 8): the exact theta-arc on which an observed Lloyd
# trajectory is realized. Ground truth = re-running deterministic Lloyd along
# the fiber and checking the trajectory matches inside the arc, differs outside.
suppressMessages(library(KmeansInference))
source("sims/unknown_var_union.R")        # fun_P0, fun_P1, fun_ts
source("sims/rfiber_exact.R")
.lloyd <- KmeansInference:::.kmeans_estimation_fixed_init

build_fiber <- function(X, cl, k1, k2) {
  P0X <- fun_P0(cl, k1, k2) %*% X; P1X <- fun_P1(cl, k1, k2) %*% X; P2X <- X - P0X - P1X
  a0 <- sqrt(sum(P0X^2)); a1 <- sqrt(sum(P1X^2)); Tt <- a0^2 + a1^2; d <- sum(cl==k1)+sum(cl==k2)-2
  list(cmat = P2X, amat = sqrt(Tt) * (P0X / a0), bmat = sqrt(Tt) * (P1X / a1),
       d = d, R_obs = d * a0^2 / a1^2)
}
xtheta <- function(fib, th) fib$cmat + sin(th) * fib$amat + cos(th) * fib$bmat
traj_at <- function(fib, th, k, init) {  # full Lloyd trajectory at fiber point th
  fr <- .lloyd(xtheta(fib, th), k, init, 10, tol_eps = 1e-6, verbose = FALSE)
  fr$cluster[seq_len(fr$iter)]
}
same_traj <- function(a, b) length(a) == length(b) && all(mapply(function(u, v) all(u == v), a, b))

npass <- 0L; nfail <- 0L
ck <- function(nm, ok, info="") { if (isTRUE(ok)) { npass<<-npass+1L; cat(sprintf("  PASS %-34s %s\n",nm,info)) } else { nfail<<-nfail+1L; cat(sprintf("  FAIL %-34s %s\n",nm,info)) } }

for (q in c(2L, 5L)) for (seed_d in c(1L, 4L, 9L)) {
  set.seed(seed_d); n_per <- 12L; k <- 3L; delta <- 4
  tc <- rep(1:k, each = n_per); mu <- rbind(c(delta,0), c(0, delta), c(-delta,0))
  X <- matrix(rnorm(n_per*k*q), n_per*k, q); X[,1:2] <- X[,1:2] + mu[tc,]
  est <- kmeans_estimation(X, k, 10, seed = 2021, verbose = FALSE)
  if (length(unique(est$final_cluster)) < k) next
  cl <- est$final_cluster; init <- est$random_init_obs; k1<-1; k2<-3
  fib <- build_fiber(X, cl, k1, k2)
  th_obs <- theta_of_r(fib$R_obs, fib$d)
  traj_obs <- traj_at(fib, th_obs, k, init)              # == observed trajectory
  arc <- path_arc(list(cmat=fib$cmat,amat=fib$amat,bmat=fib$bmat), traj_obs, init, k)
  tag <- sprintf("q=%d seed=%d", q, seed_d)
  in_arc <- nrow(arc) > 0 && any(arc[,1] <= th_obs & th_obs <= arc[,2])
  ck(paste("theta_obs in arc", tag), in_arc, sprintf("arc=%s th_obs=%.3f", paste(round(arc,3),collapse=","), th_obs))
  if (!in_arc) next
  # the component containing theta_obs
  comp <- arc[which(arc[,1] <= th_obs & th_obs <= arc[,2])[1], ]
  # INSIDE: 5 points strictly inside -> trajectory must equal observed
  ins <- seq(comp[1] + 1e-4, comp[2] - 1e-4, length.out = 5)
  ok_in <- all(vapply(ins, function(t) same_traj(traj_at(fib,t,k,init), traj_obs), logical(1)))
  ck(paste("trajectory const INSIDE arc", tag), ok_in)
  # JUST OUTSIDE: trajectory must differ (skip if component touches the boundary)
  outs <- c()
  if (comp[1] > 1e-3) outs <- c(outs, comp[1] - 1e-4)
  if (comp[2] < pi/2 - 1e-3) outs <- c(outs, comp[2] + 1e-4)
  ok_out <- length(outs) == 0 || all(vapply(outs, function(t) !same_traj(traj_at(fib,t,k,init), traj_obs), logical(1)))
  ck(paste("trajectory CHANGES outside arc", tag), ok_out)
}
cat(sprintf("\n==== %d passed, %d failed ====\n", npass, nfail))
if (nfail > 0) quit(status = 1)
