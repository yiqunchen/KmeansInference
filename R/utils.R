# ----- general purpose helper functions -----
# ----- Imported from https://github.com/lucylgao/clusterpval/blob/master/R/util.R -----
#' Takes the l2-norm of a vector.
#'
#' @keywords internal
#'
#' @param x the vector to be normed
#'
#' @return Returns the l2-norm of x.
norm_vec <- function(x) {
  ell_2_norm <- sqrt(sum(x^2))
  return(ell_2_norm)
}

#' Checks if input is an integer between a and b
#'
#' @keywords internal
#'
#' @param x input to check
#' @param a lower
#' @param b upper
#'
#' @return Returns TRUE if input is an integer between a and b, FALSE otherwise
is_integer_between_a_b <- function(x, a, b) {
  result_integer <- (x>= min(c(a, b))) && (x %% 1 == 0) && (x <= max(c(a, b)))
  return(result_integer)
}

#' Checks if two clusterings are the same up to permutation
#'
#' @keywords internal
#'
#' @param cl1 the first clustering
#' @param cl2 the second clustering
#' @param K the number of clusters
#'
#' @return Returns TRUE if they are the same, and FALSE otherwise
same_cl <- function(cl1, cl2, K) {
  tab <- table(cl1, cl2)
  same_up_to_perm <- sum(tab != 0) == K
  return(same_up_to_perm)
}

#' Checks if Ck, Ck' in C(x'(phi))
#'
#' @keywords internal
#'
#' @param cl clustering of x
#' @param cl_phi clustering of x'(phi)
#' @param k1 index of clusters involved in the test
#' @param k2 index of clusters involved in the test
#'
#' @return Returns TRUE if Ck, Ck' in C(x'(phi)), and FALSE otherwise
preserve_cl <- function(cl, cl_phi, k1, k2) {
  tab <- table(cl, cl_phi)

  k1_in <- (sum(tab[k1, ] != 0) == 1) & (sum(tab[, k1] != 0) == 1)
  k2_in <- (sum(tab[k2, ] != 0) == 1) & (sum(tab[, k2] != 0) == 1)

  return(k1_in & k2_in)
}

#' @keywords internal
#' @export
multivariate_Z_test <- function(X, cluster_vec, k1, k2, sig) {
  q <- ncol(X)
  diff_means <- colMeans(X[cluster_vec == k1, , drop=F]) -
    colMeans(X[cluster_vec == k2, , drop=F])
  stat <- norm_vec(diff_means)
  n1 <- sum(cluster_vec == k1)
  n2 <- sum(cluster_vec == k2)
  squared_norm_nu <- 1/n1 + 1/n2
  scale_factor <- squared_norm_nu*sig^2
  accurate_pchi <- pchisq(stat^2/scale_factor, df=q, log.p = TRUE,lower.tail=FALSE)
  pval <- exp(accurate_pchi) #1 - pchisq(stat^2/scale_factor, df=q)
  return(list(stat=stat, pval=pval))
}

#' Summarize the inferential result for k-means clustering
#' @param object output from running kmeans_inference
#' @param ... to be passed to methods
#' @return A data frame with summarized results
#' @export
#' @examples
#' lev1 <- 0 # mean for group 1
#' lev2 <- 3 # mean (absolute value) for group 2/3
#' sigma <- 1 # level of noise
#' nn <- 8 # grid size
#' Dmat <- genlasso::getD2d(nn, nn) # generate D matrix for the 2D fused lasso
#' ### Create the underlying signal
#' A <- matrix(lev1, ncol=nn, nrow = nn)
#' A[1:round(nn/3),1:round(nn/3)] <- 1*lev2
#' A[(nn-2):(nn),(nn-2):(nn)] <- -1*lev2
#' ### Visualize the underlying signal
#' lattice::levelplot(A)
#' set.seed(2005)
#' A.noisy <- A + rnorm(nn^2,mean=0,sd=sigma)
#' y <- c(t(A.noisy))
#' ### Now use the fusedlasso function to obtain estimated connected components after K=13
#' ### steps of the dual path algorithm
#' K = 13
#' complete_sol <- genlasso::fusedlasso(y=y,D=Dmat,maxsteps=K)
#' beta_hat <- complete_sol$beta[,K]
#' ### estimated connected components
#' estimated_CC <- complete_sol$pathobjs$i
#' estimated_CC
#' ### Run a test for a difference in means between estimated connected components 1 and 2
#' result_demo <- fusedlasso_inf(y=y, D=Dmat, c1=1, c2=2, method="K",
#' sigma=sigma, K=K, compute_ci=TRUE)
#' summary(result_demo)
summary.kmeans_inference <- function(object, ...){
  result <- data.frame(cluster_1 = object$cluster_1,
                       cluster_2 = object$cluster_2,
                       test_stats = object$test_stats,
                       p_kmeans = object$pval,
                       p_naive = object$p_naive[['pval']])
  return(result)
}
