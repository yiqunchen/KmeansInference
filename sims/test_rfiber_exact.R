#!/usr/bin/env Rscript
# Unit tests for the modular exact R-fiber solver (sims/rfiber_exact.R).
# Each module is checked against an independent ground truth.
source("sims/rfiber_exact.R")
set.seed(1)
npass <- 0L; nfail <- 0L
check <- function(name, ok, info = "") {
  if (isTRUE(ok)) { npass <<- npass + 1L; cat(sprintf("  PASS  %-38s %s\n", name, info)) }
  else            { nfail <<- nfail + 1L; cat(sprintf("  FAIL  %-38s %s\n", name, info)) }
}

cat("== Module 1: r <-> theta map ==\n")
d <- 7; th <- c(0.05, 0.3, 0.8, 1.2)
r  <- r_of_theta(th, d); th2 <- theta_of_r(r, d)
check("theta<->r round trip", max(abs(th - th2)) < 1e-12, sprintf("max err %.1e", max(abs(th - th2))))
check("r monotone in theta", all(diff(r) > 0))
check("R_obs = d tan^2(theta_obs)", abs(r_of_theta(theta_of_r(33.5, d), d) - 33.5) < 1e-9)

cat("== Module 2: polynomial eval + real roots ==\n")
check("poly_eval(1-3t+2t^2 @t=2)=3", abs(poly_eval(c(1, -3, 2), 2) - 3) < 1e-12)
rts <- poly_real_roots_in(c(0.1, -0.7, 1), 0, 1)             # (t-0.2)(t-0.5)
check("roots of t^2-0.7t+0.1 in [0,1]", length(rts) == 2 && max(abs(sort(rts) - c(0.2, 0.5))) < 1e-9,
      sprintf("got %s", paste(round(rts, 4), collapse = ",")))
# quartic (t-0.3)(t-0.7)(t^2+1): real roots 0.3,0.7 in [0,1], complex pair excluded
q4 <- c(0.21, -0.91 + 0.0, 1.0, -1.0, 1.0)  # will recompute below precisely
poly_mul <- function(a, b) { out <- numeric(length(a) + length(b) - 1)
  for (i in seq_along(a)) for (j in seq_along(b)) out[i + j - 1] <- out[i + j - 1] + a[i] * b[j]; out }
q4 <- poly_mul(poly_mul(c(-0.3, 1), c(-0.7, 1)), c(1, 0, 1))  # (t-.3)(t-.7)(t^2+1)
rq <- poly_real_roots_in(q4, 0, 1)
check("quartic real roots in [0,1] = {0.3,0.7}", length(rq) == 2 && max(abs(sort(rq) - c(0.3, 0.7))) < 1e-8,
      sprintf("got %s", paste(round(rq, 4), collapse = ",")))

cat("== Module 3+4: margin form reproduces the actual margin ==\n")
qd <- 4
rv <- function() rnorm(qd)
gj <- rv(); pj <- rv(); sj <- rv(); gk <- rv(); pk <- rv(); sk <- rv()
co <- margin_form(gj, pj, sj, gk, pk, sk)
direct <- function(al, be) sum((gj + al * pj + be * sj)^2) - sum((gk + al * pk + be * sk)^2)
ths <- runif(50, 0, pi / 2); al <- sin(ths); be <- cos(ths)
resid <- max(abs(mapply(direct, al, be) - mapply(function(a, b) eval_form(co, a, b), al, be)))
check("eval_form == direct margin", resid < 1e-12, sprintf("max resid %.1e", resid))

cat("== Module 5: N(t) identity and arc solver ==\n")
# N(t) == Q(sin th, cos th) * (1+t^2)^2  for theta = 2 atan(t)
co2 <- setNames(rnorm(6), c("A", "B", "C", "D", "E", "Fc"))
tt <- runif(50, 0, 1); thh <- 2 * atan(tt)
lhs <- poly_eval(form_to_Nt(co2), tt)
rhs <- mapply(function(a, b) eval_form(co2, a, b), sin(thh), cos(thh)) * (1 + tt^2)^2
check("N(t) == Q*(1+t^2)^2", max(abs(lhs - rhs)) < 1e-10, sprintf("max resid %.1e", max(abs(lhs - rhs))))
# arc solver vs dense theta grid (sign of Q)
ngrid <- 4000; thg <- seq(1e-6, pi / 2 - 1e-6, length.out = ngrid)
disagree <- 0
for (trial in 1:20) {
  cc <- setNames(rnorm(6), c("A", "B", "C", "D", "E", "Fc"))
  arcs <- solve_form_arcs(cc)
  in_arc <- rep(FALSE, ngrid)
  if (nrow(arcs)) for (i in seq_len(nrow(arcs))) in_arc <- in_arc | (thg >= arcs[i, 1] & thg <= arcs[i, 2])
  qle <- mapply(function(a, b) eval_form(cc, a, b), sin(thg), cos(thg)) <= 0
  # allow mismatch only immediately adjacent to a boundary
  bd <- abs(c(0, diff(in_arc))) > 0 | abs(c(diff(in_arc), 0)) > 0
  disagree <- disagree + sum(in_arc != qle & !bd)
}
check("solve_form_arcs matches grid sign", disagree == 0, sprintf("interior disagreements = %d/80000", disagree))

cat("== Module 6: interval algebra ==\n")
ii <- intersect_intervals(rbind(c(0, 2), c(3, 5)), rbind(c(1, 4)))
check("intersect", nrow(ii) == 2 && max(abs(ii - rbind(c(1, 2), c(3, 4)))) < 1e-12,
      paste(apply(round(ii, 2), 1, paste, collapse = "-"), collapse = " "))
mm <- .merge_intervals(rbind(c(0, 1), c(0.5, 2), c(3, 4)))
check("merge overlapping", nrow(mm) == 2 && max(abs(mm - rbind(c(0, 2), c(3, 4)))) < 1e-12)

cat("== Module 7: truncated-F ==\n")
# full support [0, huge] -> p == survival F(R_obs)
df1 <- 2; df2 <- 76; Robs <- 1.3
p_full <- trunc_F_p(Robs, matrix(c(0, 1e6), 1), df1, df2)
check("trunc_F over full support == survival", abs(p_full - pf(Robs, df1, df2, lower.tail = FALSE)) < 1e-6,
      sprintf("%.4f vs %.4f", p_full, pf(Robs, df1, df2, lower.tail = FALSE)))
# R_obs at left edge of the (only) interval -> p == 1
check("trunc_F, R_obs at left edge == 1", abs(trunc_F_p(2, matrix(c(2, 10), 1), df1, df2) - 1) < 1e-9)
# R_obs at right edge -> p == 0
check("trunc_F, R_obs at right edge == 0", abs(trunc_F_p(10, matrix(c(2, 10), 1), df1, df2) - 0) < 1e-9)

cat(sprintf("\n==== %d passed, %d failed ====\n", npass, nfail))
if (nfail > 0) quit(status = 1)
