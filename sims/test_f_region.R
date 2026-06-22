#!/usr/bin/env Rscript
# Unit tests for the F-REGION (unknown-variance R-fiber) solver correctness.
# Complements test_rfiber_exact.R (primitives) and test_path_arc.R (assembler) by
# checking the two pieces those don't: the truncated-F p-value math (trunc_F_p) and
# the fiber's R(theta)=d*tan^2(theta) relation -- plus the statistical pivot
# property, preservation correctness, and numerical-vs-exact agreement.
# Every check is against an INDEPENDENT ground truth (Monte Carlo, analytic
# survival ratios, re-running k-means, or a second implementation).
suppressMessages({ library(KmeansInference) })
source("sims/rfiber_exact.R")               # trunc_F_p, r_of_theta, theta_of_r
source("sims/unknown_var_union.R")          # fun_P0, fun_P1, fun_ts
source("sims/kmeans_union_unknownvar.R")        # numerical solver
source("sims/kmeans_union_unknownvar_exact.R")  # exact solver
set.seed(1)
np <- 0L; nf <- 0L
check <- function(name, ok, info="") { if (isTRUE(ok)) { np <<- np+1L; cat(sprintf("  PASS  %-44s %s\n",name,info)) }
  else { nf <<- nf+1L; cat(sprintf("  FAIL  %-44s %s\n",name,info)) } }

## ===================================================================
cat("== A. trunc_F_p vs Monte-Carlo (the truncated-F p-value formula) ==\n")
df1 <- 5L; df2 <- 40L
Fs <- rf(4e6, df1, df2)
mc <- function(R, ivl) { inS <- rep(FALSE, length(Fs))
  for (i in seq_len(nrow(ivl))) inS <- inS | (Fs >= ivl[i,1] & Fs <= ivl[i,2])
  sum(Fs >= R & inS) / sum(inS) }
cases <- list(
  list(nm="single [0.5,Inf)",   R=2.0, ivl=matrix(c(0.5, Inf), 1)),
  list(nm="bounded [0.8,3.0]",   R=1.5, ivl=matrix(c(0.8, 3.0), 1)),
  list(nm="two arcs",            R=2.0, ivl=matrix(c(0.4,1.0, 1.6,3.5), 2, byrow=TRUE)),
  list(nm="R in 2nd arc",        R=2.2, ivl=matrix(c(0.4,1.0, 1.6,3.5), 2, byrow=TRUE)),
  list(nm="three arcs",          R=1.2, ivl=matrix(c(0.3,0.7, 1.0,1.5, 2.0,4.0), 3, byrow=TRUE)))
for (c. in cases) { got <- trunc_F_p(c.$R, c.$ivl, df1, df2); want <- mc(c.$R, c.$ivl)
  check(sprintf("MC: %s", c.$nm), abs(got-want) < 0.01, sprintf("got=%.4f mc=%.4f", got, want)) }

cat("== B. trunc_F_p analytic survival ratio + tail/underflow + edges ==\n")
# single [a,Inf): p = S(R)/S(a) exactly (S=upper tail)
a <- 0.8; R <- 2.5; want <- pf(R,df1,df2,lower.tail=FALSE)/pf(a,df1,df2,lower.tail=FALSE)
check("analytic S(R)/S(a)", abs(trunc_F_p(R, matrix(c(a,Inf),1), df1, df2) - want) < 1e-10,
      sprintf("=%.6f", want))
# extreme tail (where naive CDF-differencing underflows): big df, huge R
dft <- 10L; dd <- 200L; Rbig <- 60; abig <- 50
wbig <- pf(Rbig,dft,dd,lower.tail=FALSE)/pf(abig,dft,dd,lower.tail=FALSE)
gbig <- trunc_F_p(Rbig, matrix(c(abig,Inf),1), dft, dd)
check("extreme-tail not NA + matches survival", is.finite(gbig) && abs(gbig-wbig) < 1e-9,
      sprintf("got=%.3e want=%.3e", gbig, wbig))
check("R below all arcs -> p=1",  abs(trunc_F_p(0.2, matrix(c(0.5,2.0),1), df1,df2) - 1) < 1e-12)
check("R above all arcs -> p=0",  abs(trunc_F_p(5.0, matrix(c(0.5,2.0),1), df1,df2) - 0) < 1e-12)

cat("== C. truncated-F is a valid PIVOT: p ~ Uniform under fixed truncation ==\n")
S <- matrix(c(0.3,0.9, 1.4,2.2, 3.0,Inf), 3, byrow=TRUE)  # arbitrary fixed set
samp <- rf(3e5, df1, df2); inS <- rep(FALSE,length(samp))
for (i in seq_len(nrow(S))) inS <- inS | (samp>=S[i,1] & samp<=S[i,2])
ps <- vapply(samp[inS], function(f) trunc_F_p(f, S, df1, df2), numeric(1))
ksp <- ks.test(ps, "punif")$p.value
check("pivot uniform (KS p>0.05)", ksp > 0.05, sprintf("KSp=%.3f  rej@.05=%.3f (n=%d)", ksp, mean(ps<=0.05), length(ps)))

## ===================================================================
cat("== D. fiber consistency: R(theta)=d*tan^2(theta) and x(theta_obs)=X ==\n")
mk <- function(seed, q=2, nper=15, delta=5) { set.seed(seed)
  mu <- rbind(c(0,0), c(delta,0), c(delta/2, delta*sqrt(3)/2)); mu <- cbind(mu, matrix(0,3,q-2))
  matrix(rnorm(3*nper*q), 3*nper, q) + mu[rep(1:3,each=nper), ] }
X <- mk(11); est <- kmeans_estimation(X, 3, 10, 2021, verbose=FALSE); cl <- est$final_cluster
k1<-1; k2<-2; n<-nrow(X)
P0 <- fun_P0(cl,k1,k2); P1 <- fun_P1(cl,k1,k2)
P0X <- P0%*%X; P1X <- P1%*%X; P2X <- X - P0X - P1X
Tt <- sum(P0X^2)+sum(P1X^2); m <- sum(cl==k1)+sum(cl==k2); d <- m-2
U0 <- P0X/sqrt(sum(P0X^2)); U1 <- P1X/sqrt(sum(P1X^2))
xth <- function(th) P2X + sqrt(Tt)*(sin(th)*U0 + cos(th)*U1)
R_obs <- fun_ts(X, cl, k1, k2); th_obs <- atan(sqrt(R_obs/d))
check("x(theta_obs) == X", max(abs(xth(th_obs) - X)) < 1e-8, sprintf("maxdiff=%.1e", max(abs(xth(th_obs)-X))))
check("fun_ts(X) == d*tan^2(theta_obs)", abs(R_obs - d*tan(th_obs)^2) < 1e-8)
for (th in c(0.25, 0.55, 0.85, 1.15)) {
  rr <- fun_ts(xth(th), cl, k1, k2)
  check(sprintf("R(theta=%.2f) == d*tan^2", th), abs(rr - r_of_theta(th, d)) < 1e-6,
        sprintf("got=%.4f want=%.4f", rr, r_of_theta(th, d))) }

cat("== E. solver sanity: theta_obs preserved, R_obs in the set ==\n")
num <- kmeans_union_unknownvar(X, 3, k1, k2, seed=2021, n_theta=600)
check("R_obs in preserved set", any(num$intervals_r[,1] <= num$R_obs & num$R_obs <= num$intervals_r[,2]),
      sprintf("R_obs=%.3f n_int=%d", num$R_obs, num$n_intervals))
check("R_obs(solver) == fun_ts(X)", abs(num$R_obs - R_obs) < 1e-6)

cat("== F. preservation correctness: inside arcs preserved, gaps NOT ==\n")
.lloyd <- KmeansInference:::.kmeans_estimation_fixed_init
.presrv <- KmeansInference:::.kmeans_preserves_observed_clusters
init <- est$random_init_obs
preserved <- function(th) { fr <- tryCatch(.lloyd(xth(th),3,init,10,tol_eps=1e-6,verbose=FALSE), error=function(e) NULL)
  if (is.null(fr)) FALSE else isTRUE(.presrv(cl, fr$final_cluster, k1, k2)) }
ivlT <- cbind(theta_of_r(num$intervals_r[,1], d), theta_of_r(num$intervals_r[,2], d))
inside_ok <- all(vapply(seq_len(nrow(ivlT)), function(i) preserved(mean(ivlT[i,])), logical(1)))
check("interior of every arc is preserved (delta=5)", inside_ok)
# find a MULTI-ARC dataset (global null splits noise -> several preserved arcs) so
# the gap-not-preserved direction is actually exercised, not skipped.
multi <- NULL
for (s in 1:40) { set.seed(2000+s); Xm <- matrix(rnorm(45*2), 45, 2)
  em <- tryCatch(kmeans_estimation(Xm,3,10,2021,verbose=FALSE), error=function(e) NULL); if (is.null(em)) next
  nm <- tryCatch(kmeans_union_unknownvar(Xm,3,1,2,seed=2021,n_theta=800), error=function(e) NULL)
  if (!is.null(nm) && nm$n_intervals >= 2) { multi <- list(X=Xm, est=em, num=nm); break } }
if (is.null(multi)) check("multi-arc gap test", FALSE, "(no multi-arc dataset found)") else {
  Xm<-multi$X; em<-multi$est; clm<-em$final_cluster; initm<-em$random_init_obs
  P0m<-fun_P0(clm,1,2); P1m<-fun_P1(clm,1,2); P0Xm<-P0m%*%Xm; P1Xm<-P1m%*%Xm; P2Xm<-Xm-P0Xm-P1Xm
  Ttm<-sum(P0Xm^2)+sum(P1Xm^2); dm<-sum(clm==1)+sum(clm==2)-2
  U0m<-P0Xm/sqrt(sum(P0Xm^2)); U1m<-P1Xm/sqrt(sum(P1Xm^2))
  xthm<-function(th) P2Xm + sqrt(Ttm)*(sin(th)*U0m+cos(th)*U1m)
  presm<-function(th){ fr<-tryCatch(.lloyd(xthm(th),3,initm,10,tol_eps=1e-6,verbose=FALSE),error=function(e)NULL)
    if(is.null(fr))FALSE else isTRUE(.presrv(clm,fr$final_cluster,1,2)) }
  ivTm<-cbind(theta_of_r(multi$num$intervals_r[,1],dm), theta_of_r(multi$num$intervals_r[,2],dm))
  ivTm<-ivTm[order(ivTm[,1]),,drop=FALSE]
  ins<-all(vapply(seq_len(nrow(ivTm)), function(i) presm(mean(ivTm[i,])), logical(1)))
  gaps<-vapply(seq_len(nrow(ivTm)-1), function(i){ g<-(ivTm[i,2]+ivTm[i+1,1])/2
    if(ivTm[i+1,1]-ivTm[i,2] < 1e-3) NA else presm(g) }, logical(1))
  gaps<-gaps[!is.na(gaps)]
  check(sprintf("multi-arc (%d arcs): all interiors preserved", nrow(ivTm)), ins)
  check("multi-arc: gap points NOT preserved", length(gaps)>0 && !any(gaps),
        sprintf("%d gaps tested", length(gaps))) }

cat("== G. numerical solver vs EXACT solver agree ==\n")
ex <- tryCatch(kmeans_union_unknownvar_exact(X, 3, k1, k2, seed=2021), error=function(e) NULL)
if (is.null(ex)) check("exact solver runs", FALSE) else {
  check("exact R_obs == numerical R_obs", abs(ex$R_obs - num$R_obs) < 1e-6)
  check("exact p_union ~ numerical p_union", abs(ex$p_union - num$p_union) < 0.02,
        sprintf("exact=%.4f num=%.4f", ex$p_union, num$p_union)) }

cat(sprintf("\n==== F-region solver: %d passed, %d failed ====\n", np, nf))
if (nf > 0) quit(status = 1)
