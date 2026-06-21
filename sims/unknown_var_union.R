#!/usr/bin/env Rscript
# ===========================================================================
# UNKNOWN-VARIANCE union / path k-means selective test (works for ANY K).
#
# Background (Yun & Barber 2023, arXiv:2301.12999; repo functions.R):
#   The unknown-variance machinery is CLUSTERING-AGNOSTIC. Given a generic
#   KNOWN-VARIANCE truncation set S (in the phi = ||x^T nu|| scale), we
#   studentize it, S' = S / (||P1 X|| * ||v||), map the resulting F-region to
#   chi^2 via SFA (Li & Martin 2002), and feed it to the SAME truncated-chi^2
#   ratio approximation (TChisqRatioApprox) that the known-variance test uses.
#   The reference distribution is F_{q, (m-2)q} with
#       m = |C_{k1}| + |C_{k2}|     (the two tested clusters ONLY).
#
#   CRITICAL FIX vs the closed-form fun_pval_sig_unknown in their repo, which
#   hardcodes n = length(cl): that only equals m when K = 2. For K > 2 the two
#   tested clusters are a strict subset, so we MUST use m = |C_{k1}| + |C_{k2}|
#   everywhere (numerator df = q, denominator df = (m-2)*q). This module does
#   that, so it is correct for any K >= 2.
#
# NULL ASSUMPTION -- H0' (POINTWISE):
#   The studentization denominator ||P1 X|| is central (its expectation is the
#   noise scale, not contaminated by signal) only under the STRONGER null
#       H0': for all i in C_{k1} u C_{k2}, the mean mu_i equals the common
#            cluster mean of its tested cluster (pointwise equality of means
#            WITHIN each of the two tested clusters),
#   which is stronger than the usual H0 (equal cluster means: mu_{C_k1} =
#   mu_{C_k2}). This is the price of not knowing sigma. The effect of an H0'
#   violation depends on its DIRECTION (verified by adversarial simulation):
#     * Broad-null-but-not-pointwise (equal cluster means, e.g. canceling
#       +d/-d within-cluster sub-populations): the heterogeneity only ADDS
#       energy to the within-cluster residual ||P1 X||, shrinking F -> the test
#       is CONSERVATIVE (Type-I-safe), not anti-conservative.
#     * Anti-conservative ONLY when within-cluster heterogeneity DEFLATES
#       ||P1 X|| while the numerator carries genuine between-cluster signal --
#       i.e. the cluster means actually differ (H0 itself is false, a non-null
#       regime). In the realistic k-means pipeline this is hard to trigger
#       because k-means tends to spend an extra cluster resolving within-cluster
#       heterogeneity, restoring H0' at fit time.
# ===========================================================================

## --- self-contained helpers (copied from the K=2 prototype) ----------------
# Contrast vector v selecting the difference in means between clusters k1, k2.
fun_v  <- function(cl, k1=1, k2=2) (cl==k1)/sum(cl==k1) - (cl==k2)/sum(cl==k2)

# Rank-1 projection onto span(v) (numerator / signal direction).
fun_P0 <- function(cl, k1=1, k2=2){ v <- fun_v(cl,k1,k2); tcrossprod(v)/sum(v^2) }

# Within-tested-cluster centering projection (denominator / noise direction):
# it removes each tested cluster's mean, leaving the within-cluster residuals.
fun_P1 <- function(cl, k1=1, k2=2){
  n1<-sum(cl==k1); n2<-sum(cl==k2); i1<-as.numeric(cl==k1); i2<-as.numeric(cl==k2)
  diag(i1)-tcrossprod(i1)/n1 + diag(i2)-tcrossprod(i2)/n2
}

# Observed F statistic: (||P0 X||^2 / q) / (||P1 X||^2 / ((m-2)q)), m=|C_k1|+|C_k2|.
fun_ts <- function(X, cl, k1=1, k2=2){
  m<-sum(cl==k1)+sum(cl==k2); q<-ncol(X)
  P0X<-fun_P0(cl,k1,k2)%*%X; P1X<-fun_P1(cl,k1,k2)%*%X
  (sum(P0X^2)/q)/(sum(P1X^2)/((m-2)*q))
}

# ||P1 X||_F : the studentization scale.
fun_P1X_norm <- function(X, cl, k1=1, k2=2) sqrt(sum((fun_P1(cl,k1,k2)%*%X)^2))

# SFA F->chi^2 map (Li & Martin 2002), elementwise on a 2-col interval matrix.
# n1, n2 are the numerator / denominator degrees of freedom of the F.
fun_F_to_chi2 <- function(x, n1, n2){
  out <- x
  fin <- is.finite(x)
  t1 <- 2*n2 + (n1*x)/3 + (n1-2); t2 <- 2*n2 + (4*n1*x)/3
  out[fin] <- ((t1/t2)*n1*x)[fin]; out[!fin] <- Inf; out
}

## --- generic unknown-variance p-value from a known-variance set S ----------
#' Unknown-variance truncated-F p-value from a GENERIC known-variance set.
#'
#' @param X  n x q data matrix.
#' @param cl integer cluster labels (length n).
#' @param S  an intervals::Intervals over phi (a known-variance truncation set,
#'           e.g. fit$interval_union or fit$interval_path).
#' @param k1,k2 the two tested cluster labels.
#' @return a scalar p-value in [0,1], or NA if S is empty.
#'
#' Uses m = |C_k1| + |C_k2| (NOT length(cl)); df_denom = (m-2)*q. Correct
#' for any K >= 2. Valid under H0' (pointwise equality of means; see header).
pval_union_unknown <- function(X, cl, S, k1=1, k2=2){
  if(is.null(S) || nrow(S)==0) return(NA_real_)
  m <- sum(cl==k1) + sum(cl==k2)         # tested clusters only -- the K>2 fix
  q <- ncol(X)
  P1Xn <- fun_P1X_norm(X, cl, k1, k2)
  vn   <- sqrt(sum(fun_v(cl, k1, k2)^2))
  Sm   <- as.matrix(S) / (P1Xn * vn)     # studentize: S' = S / (||P1 X|| ||v||)
  den  <- (m - 2) * Sm^2                 # known-var phi set -> F-region (>= 0)
  ts   <- fun_ts(X, cl, k1, k2)          # observed F statistic
  num  <- suppressWarnings(intervals::interval_intersection(
            intervals::Intervals(c(ts, Inf)),
            intervals::Intervals(den)))
  if(nrow(num)==0) return(0)
  den_c <- intervals::Intervals(fun_F_to_chi2(den,            q, (m-2)*q))
  num_c <- intervals::Intervals(fun_F_to_chi2(as.matrix(num), q, (m-2)*q))
  KmeansInference:::TChisqRatioApprox(q, num_c, den_c)
}

## --- convenience wrapper around a kmeans_inference_union fit ----------------
#' Unknown-variance PATH p-value from a fitted union object.
#'
#' Only the PATH unknown-variance test is returned, because it is the only valid
#' one: feeding the UNION set into the studentized-F remap is ANTI-CONSERVATIVE
#' (removed). The studentized-F pivot is valid only conditional on a SINGLE Lloyd
#' path; the union conditions on less and the F-pivot does not survive that
#' weaker conditioning (verified EXP6: exact-F union over-rejects 0.081 and
#' scales with union width 0.010->0.170, while the exact-chi union on the SAME
#' set is calibrated 0.037; EXP7: F-path 0.047/0.067 at q=2/10). The PATH test
#' reproduces Yun-Barber and is calibrated under H0'.
#'
#' To get a VALID more-powerful union under unknown variance you must explore in
#' the F-pivot's OWN coordinate R (perturb along x'(r) that co-varies P0X,P1X at
#' fixed total energy), NOT remap the phi-union set -- see sims/fix_xr_validation.R.
#'
#' @param fit    output of KmeansInference::kmeans_inference_union.
#' @param X      the data matrix that produced `fit`.
#' @param k1,k2  the two tested cluster labels (default to fit$cluster_1/_2).
#' @return list(p_path_unknown = [valid Yun-Barber path test]).
unknown_var_from_fit <- function(fit, X, k1=NULL, k2=NULL){
  if(is.null(k1)) k1 <- fit$cluster_1
  if(is.null(k2)) k2 <- fit$cluster_2
  cl <- fit$final_cluster
  list(p_path_unknown = pval_union_unknown(X, cl, fit$interval_path, k1, k2))
}
