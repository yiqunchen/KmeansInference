#' Perform k-means clustering on a data matrix.
#'
#' @param X Numeric matrix; \eqn{n} by \eqn{q} matrix of observed data
#' @param k Integer; the number of clusters for k-means clustering
#' @param iter.max Positive integer; 	the maximum number of iterations allowed in k-means clustering (Lloyd's) algorithm.
#' Default to \code{10}.
#' @param seed Random seed for the initialization in k-means clustering algorithm.
#'
#' @details
#' The data given by X are clustered by k-means clustering,
#' which aims to partition the points into k groups such that the sum of squares from points
#' to the assigned cluster centers is minimized. In other words, k-means clustering solves
#' the following optimization problem
#' \deqn{ \sum_{k=1}^K \sum_{i \in \mathcal{C}_k} \left\Vert x_i -  \frac{\sum_{i \in \mathcal{C}_k} x_i}{|\mathcal{C}_k|}
#'  \right\Vert_2^2 , }
#'  subject the constraint that \eqn{\mathcal{C}_1,..., {\mathcal{C}_K}} forms a partition of the integers \eqn{1,..., n}.
#' The algorithm from Lloyd (1957) (also proposed in MacQueen (1967)) is used to produce a solution.
#'
#' This function is a re-implementation of the kmeans function in base R (i.e., the stats package) that
#' stores all the intermediate clustering assignments as well (see Section 3 of our manuscript for details).
#' Ouputs from these two functions agree on their estimated clusters, as well as their ordering.
#'
#' N.B.: the kmeans function in base R was implemented in Fortran and C, while our implementation is entirely in R.
#' As a result, there might be corner cases where these two functions disagree.
#' @return Returns a list with the following elements:
#' \itemize{
#' \item \code{final_cluster} Estimated clusters via k-means clustering
#' \item \code{centers} A matrix of the cluster centroids.
#' \item \code{objective} The objective function at the final iteration of k-means algorithm.
#' }
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
#' @references
#' Lloyd, S. P. (1957, 1982). Least squares quantization in PCM. Technical Note, Bell Laboratories.
#' Published in 1982 in IEEE Transactions on Information Theory, 28, 128–137.
#'
#' MacQueen, J. (1967). Some methods for classification and analysis of multivariate observations.
#' In Proceedings of the Fifth Berkeley Symposium on Mathematical Statistics and Probability,
#' pp. 281–297. Berkeley, CA: University of California Press.
#' @export
kmeans_estimation <- function(X, k, iter.max = 10, seed = 1234,
                              tol_eps = 1e-4, verbose=TRUE){

  set.seed(seed)
  if(!is.matrix(X)) stop("X should be a matrix")
  if(k>=nrow(X)){
    stop("Cannot have more clusters than observations")
  }
  iter_T <- 0
  n <- dim(X)[1]
  p <- dim(X)[2]
  cluster_assign_list <- vector("list", length = iter.max)
  centroid_list <- vector("list", length = iter.max)
  objective_value <- vector("list", length = iter.max)
  # first set of centroids
  initial_sample <- sample(c(1:n),k,replace=F)
  current_centroid <- X[initial_sample,]
  # first set of assignments
  distance_matrix <- rdist::cdist(current_centroid,X)
  current_cluster <- apply(distance_matrix,2,which.min)
  iter_T <- iter_T+1
  centroid_list[[iter_T]] <- current_centroid
  cluster_assign_list[[iter_T]] <- current_cluster
  curr_objective_value <- sum(apply(distance_matrix,2,min)^2)
  objective_diff <- 10000 #curr_objective_value
  objective_value[[iter_T]] <- curr_objective_value
  same_cluster <- FALSE
  while((iter_T<=iter.max)&(!same_cluster)){
    # update centroids
    for (current_k in c(1:k)){
      X_current <- X[(current_cluster==current_k), ,drop=F]
      new_centroid_k <- .colMeans(X_current, dim(X_current)[1], dim(X_current)[2])
      current_centroid[current_k,] <- new_centroid_k  # 1 by q
    } # current_centroid is k by q
    # update assignments
    distance_matrix <- rdist::cdist(current_centroid,X)
    current_cluster <- apply(distance_matrix,2,which.min)
    # add iteration and store relevant information
    iter_T <- iter_T+1
    centroid_list[[iter_T]] <- current_centroid
    cluster_assign_list[[iter_T]] <- current_cluster
    same_cluster <- all(current_cluster==cluster_assign_list[[iter_T-1]])
    # update objective diff
    new_objective_value <- sum(apply(distance_matrix,2,min)^2)
    objective_diff <- abs(curr_objective_value-(new_objective_value))/curr_objective_value
    curr_objective_value <- new_objective_value
    # store objextive as well
    objective_value[[iter_T]] <- curr_objective_value
  }

  result_list <- list("cluster" = cluster_assign_list, "centers" = centroid_list,
                      "objective" = objective_value, "iter" = iter_T,
                      "final_cluster" = cluster_assign_list[[iter_T]],
                      "random_init_obs" = initial_sample)
  return(result_list)
}


# ----- main function to test equality of the means of two estimated clusters via k-means clustering -----
#' Test for a difference in means between clusters of observations
#' identified via k-means clustering.
#'
#' This function tests the null hypothesis of no difference in means between
#' output by k-means clustering. The clusters are numbered as per the results of
#' the \code{kmeans_estimation} function in the \code{KmeansInference} package.
#' @param X Numeric matrix; \eqn{n} by \eqn{q} matrix of observed data
#' @param k Integer; the number of clusters for k-means clustering
#' @param cluster_1,cluster_2 Two different integers in {1,...,k}; two estimated clusters to test, as indexed by the results of
#' \code{kmeans_estimation}.
#' @param iso Boolean. If TRUE, an isotropic covariance matrix model is used.
#' @param sig Numeric; noise standard deviation for the observed data, a non-negative number;
#' relevant if \code{iso}=TRUE. If it's not given as input, a median-based estimator will be by default (see Section 4.2 of our manuscript).
#' @param SigInv Numeric matrix; if \code{iso} is FALSE, *required* \eqn{q} by \eqn{q} matrix specifying \eqn{\Sigma^{-1}}.
#' @param iter.max Positive integer; 	the maximum number of iterations allowed in k-means clustering algorithm. Default to \code{10}.
#' @param seed Random seed for the initialization in k-means clustering algorithm.
#' @param tol_eps A small number specifying the convergence criterion for k-means clustering,
#' default to \code{1e-6}.
#'
#' @return Returns a list with the following elements:
#' \itemize{
#' \item \code{pval} the selective p-value \eqn{p_{selective}} in Chen and Witten (2022+)
#' \item \code{final_interval} the conditioning set of Chen and Witten (2022+), stored as the \code{Intervals} class
#' \item \code{test_stats} test statistic: the difference in the empirical means of two estimated clusters
#' \item \code{final_cluster} Estimated clusters via k-means clustering
#' }
#'
#' @export
#'
#' @details
#' Consider the generative model \eqn{X \sim MN(\mu,I_n,\sigma^2 I_q)}, k-means clustering
#' solves the following optimization problem
#' \deqn{ \sum_{k=1}^K \sum_{i \in \mathcal{C}_k} \left\Vert x_i -  \frac{\sum_{i \in \mathcal{C}_k} x_i}{|\mathcal{C}_k|}
#'  \right\Vert_2^2 , }
#'  where \eqn{\mathcal{C}_1,..., {\mathcal{C}_K}} forms a partition of the integers \eqn{1,..., n}, and can be regarded as
#'  the estimated clusters of the original observations. In practice, solutions to the optimization problem is
#'  often obtained using iterative algorithms, e.g., the Lloyd's algorithm.
#' Now suppose we want to test whether the means of two estimated clusters \code{cluster_1} and \code{cluster_2}
#' are equal; or equivalently, the null hypothesis of the form \eqn{H_{0}:  \mu^T \nu = 0_q} versus
#' \eqn{H_{1}:   \mu^T \nu \neq 0_q} for suitably chosen \eqn{\nu} and all-zero vectors \eqn{0_q}.
#'
#' This function computes the following p-value:
#' \deqn{P \left( \left\Vert X^T \nu \right\Vert_2 \ge \left\Vert x^T \nu  \right\Vert_2 \; | \;
#'   \bigcap_{t=1}^{T}\bigcap_{i=1}^{n} \left\{ c_i^{(t)} \left( X \right) =
#'  c_i^{(t)}\left( x \right) \right\},  \Pi_\nu^\perp Y  =  \Pi_\nu^\perp y \right),}
#' where \eqn{c_i^{(t)}} is the is the cluster assigned to the \eqn{i}th observation at the \eqn{t}th iteration of
#' the Lloyd's algorithm, and \eqn{\Pi_\nu^\perp} is the orthogonal projection to the orthogonal complement of \eqn{\nu}.
#' In particular, the test based on this p-value controls the selective Type I error and has substantial power.
#' Readers can refer to the Sections 2 and 4 in Chen and Witten (2022+) for more details.
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
#' @references
#' Chen YT, Witten DM. (2022+) Selective inference for k-means clustering. arXiv preprint.
#' https://arxiv.org/abs/xxxx.xxxxx.
#' Lloyd, S. P. (1957, 1982). Least squares quantization in PCM. Technical Note, Bell Laboratories.
#' Published in 1982 in IEEE Transactions on Information Theory, 28, 128–137.
#'
kmeans_inference <- structure(function(X, k, cluster_1, cluster_2,
                                       iso=TRUE, sig=NULL, SigInv=NULL,
                             iter.max = 10, seed = 1234,
                             tol_eps = 1e-6, verbose=TRUE){

  set.seed(seed)
  if(!is.matrix(X)) stop("X should be a matrix")
  if(k>=nrow(X)){
    stop("Cannot have more clusters than observations")
  }
  if(is.null(sig)&is.null(SigInv)){
    stop("At least one of variance and covariance matrix must be specified!")
  }
  if((!is.null(sig))&(!is.null(SigInv))){
    stop("Only one of variance and covariance matrix can be specified!")
  }
  if ((iso)&(is.null(sig))){
    cat("Specifying sig is needed when iso=TRUE!\n")
    cat("variance  not specified, using a robust estimator by default!\n")
    estimate_MED <- function(X){
      for (j in c(1:ncol(X))){
        X[,j] <- X[,j]-median(X[,j])}
      sigma_hat <- sqrt(median(X^2)/qchisq(1/2,df=1))
      return(sigma_hat)
    }
    sig <- estimate_MED(X)
  }
  if (!(iso)&(is.null(SigInv))){
    stop("Specifying SigInv is needed when iso=FALSE!\n")
  }
  if((min(cluster_1,cluster_2)<1)|(max(cluster_1,cluster_2)>k)){
    stop("Cluster numbers must be between 1 and k!")
  }
  n <- dim(X)[1]
  p <- dim(X)[2]
  # get the list of all assigned clusters first
  estimated_k_means <- kmeans_estimation(X, k, iter.max, seed, tol_eps, verbose)
  # check if we get the desired number of clusters:
  if(length(unique(estimated_k_means$final_cluster))<k){
    stop("k-means clustering did not return the desired number of clusters! Try a different seed?")
  }
  estimated_final_cluster <- estimated_k_means$cluster[[estimated_k_means$iter]]
  all_T_clusters <- do.call(rbind, estimated_k_means$cluster)
  all_T_centroids <- estimated_k_means$centers
  T_length <- nrow(all_T_clusters)
  # construct contrast vector
  v_vec <- rep(0, times=nrow(X))
  v_vec[estimated_final_cluster == cluster_1] = 1/(sum(estimated_final_cluster == cluster_1))
  v_vec[estimated_final_cluster == cluster_2] = -1/(sum(estimated_final_cluster == cluster_2))

  n1 <- sum(estimated_final_cluster == cluster_1)
  n2 <- sum(estimated_final_cluster == cluster_2)
  squared_norm_nu <- 1/n1 + 1/n2
  v_norm <- sqrt(squared_norm_nu) # recycle this computed value
  # compute XTv
  diff_means <- colMeans(X[estimated_final_cluster == cluster_1, ,drop=FALSE]) -
    colMeans(X[estimated_final_cluster == cluster_2, , drop=FALSE])
  # compute
  XTv <- diff_means
  XTv_norm <- norm_vec(diff_means)
  dir_XTv <- XTv/XTv_norm

  p_naive <- NULL
  # compute test_stats in the isotropic case
  if(!is.null(sig)){
    test_stats <- XTv_norm
    scale_factor <- squared_norm_nu*sig^2
    gestat <- intervals::Intervals(c(test_stats^2/scale_factor, Inf))

    # compute S
    final_interval_chisq <- kmeans_compute_S_iso(X, estimated_k_means, all_T_clusters,
                                                 all_T_centroids,
                                           n, XTv, XTv_norm,
                                           dir_XTv, v_vec,
                                           v_norm, T_length, k)

    p_naive <- multivariate_Z_test(X, estimated_final_cluster, cluster_1, cluster_2, sig)
  }

  # compute test_stats in the general cov case
  if(!is.null(SigInv)){
    test_stats <- sqrt(as.numeric(t(diff_means)%*%SigInv%*%diff_means))
    scale_factor <- squared_norm_nu
    gestat <- intervals::Intervals(c(test_stats^2/scale_factor, Inf))
    # compute S
    final_interval_chisq <- kmeans_compute_S_genCov(X, estimated_k_means, all_T_clusters,
                                                    all_T_centroids,
                                           n, XTv, XTv_norm,
                                           dir_XTv, v_vec,
                                           v_norm, T_length, test_stats, k)

  }


  # improve numerical stability
  final_interval_chisq <- intervals::interval_union(final_interval_chisq,
                           intervals::Intervals_full(c(test_stats-(1e-09),test_stats+(1e-09))),
                           check_valid=FALSE)


  denom <- final_interval_chisq^2/scale_factor

  #cat("gestat",gestat,"\n")
  #cat("denom",denom,"\n")

  numer <- suppressWarnings(intervals::interval_intersection(gestat, denom))
  pval <- TChisqRatioApprox(p, numer, denom)


  result_list <- list("final_interval"=final_interval_chisq,
                      "final_cluster" = estimated_final_cluster,
                      "test_stat"=test_stats,
                      "cluster_1" = cluster_1,
                      "cluster_2" = cluster_2,
                      "sig" = sig, "SigInv" = SigInv,
                      "scale_factor" = scale_factor,
                      "p_naive" = p_naive,
                       "call" = match.call(),
                      "pval" = pval)
  class(result_list) <- "kmeans_inference"
  return(result_list)

})







