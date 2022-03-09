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
#' library(KmeansInference)
#' library(ggplot2)
#' set.seed(2022)
#' n <- 150
#' true_clusters <- c(rep(1, 50), rep(2, 50), rep(3, 50))
#' delta <- 10
#' q <- 2
#' mu <- rbind(c(delta/2,rep(0,q-1)),
#' c(rep(0,q-1), sqrt(3)*delta/2),
#' c(-delta/2,rep(0,q-1)) )
#' sig <- 1
#' # Generate a matrix normal sample
#' X <- matrix(rnorm(n*q, sd=sig), n, q) + mu[true_clusters, ]
#' # Visualize the data
#' ggplot(data.frame(X), aes(x=X1, y=X2)) +
#' geom_point(cex=2) + xlab("Feature 1") + ylab("Feature 2") +
#'  theme_classic(base_size=18) + theme(legend.position="none") +
#'  scale_colour_manual(values=c("dodgerblue3", "rosybrown", "orange")) +
#'  theme(legend.title = element_blank(),
#'  plot.title = element_text(hjust = 0.5))
#'  k <- 3
#'  # Run k-means clustering with K=3
#'  estimated_clusters <- kmeans_estimation(X, k,iter.max = 20,seed = 2021)$final_cluster
#'  table(true_clusters,estimated_clusters)
#'  # Visualize the clusters
#'  ggplot(data.frame(X), aes(x=X1, y=X2, col=as.factor(estimated_clusters))) +
#'  geom_point(cex=2) + xlab("Feature 1") + ylab("Feature 2") +
#'  theme_classic(base_size=18) + theme(legend.position="none") +
#'  scale_colour_manual(values=c("dodgerblue3", "rosybrown", "orange")) +
#'  theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5))
#'  ### Run a test for a difference in means between estimated clusters 1 and 3
#'  cluster_1 <- 1
#'  cluster_2 <- 3
#'  cl_1_2_inference_demo <- kmeans_inference(X, k=3, cluster_1, cluster_2,
#'  sig=sig, iter.max = 20, seed = 2021)
#'  summary(cl_1_2_inference_demo)
summary.kmeans_inference <- function(object, ...){
  result <- data.frame(cluster_1 = object$cluster_1,
                       cluster_2 = object$cluster_2,
                       test_stats = object$test_stats,
                       p_selective = object$pval,
                       p_naive = object$p_naive[['pval']])
  return(result)
}
