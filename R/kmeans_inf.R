#' k means estimation
#'
#' @export
kmeans_estimation <- function(X, k, iter.max = 10, nstart = 1, seed = 1234,
                              tol_eps = 1e-2, verbose=TRUE){

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


# ----- main function to test equality of the means of two estimated connected components -----
#' Testing for a difference in means between clusters of observations
#' identified via k-means clustering.
#'
#' This functions tests the null hypothesis of no difference in means between
#' two estimated clusters \code{cluster_1} and \code{cluster_2} of the output of the
#' k means clustering solution obtained via the Lloyd's algorithm.
#' The ordering are numbered as per the results of the \code{fusedlasso}
#' function in the \code{genlasso} package.
#'
#'(X, k, cluster_1, cluster_2, sig=NULL, SigInv=NULL,
#iter.max = 10, nstart = 1,
#seed = 1234)
#' Input:
#' @param y Numeric vector; \eqn{n} dimensional observed data
#' @param D Numeric matrix; \eqn{m} by \eqn{n} penalty matrix, i.e.,
#' the oriented incidence matrix over the underlying graph
#' @param c1,c2 Integers selecting the two connected components to test, as indexed by the results of
#' \code{genlasso::fusedlasso}.
#' @param method One of "K" or "CC", which indicates which conditioning set to use
#' @param sigma Numeric; noise standard deviation for the observed data, a non-negative number.
#' @param K Integer; number of steps to run the dual-path algorithm.
#' It must be specified if method=="K".
#' @param L Integer; the targeted number of connected components.
#' It must be specified if method=="CC".
#' @param early_stop Numeric; specify when the truncation set computation
#' should be terminated. The default is NULL, which indicates infinity.
#' @param compute_ci Logical; the default is False. Specifying whether confidence intervals for \eqn{\nu^{T}\beta}, the
#' difference in means between the two estimated connected components, should be computed.
#' @param alpha_level Numeric; parameter for the 1-\code{alpha_level} confidence interval, defeault to 0.05
#' @return Returns a list with elements:
#' \itemize{
#' \item \code{pval} the p-value in Chen et al. (2021+)
#' \item \code{truncation_set} the conditioning set of Chen et al. (2021+) stored as \code{Intervals} class
#' \item \code{test_stats} test statistics: the difference in means of two connected components
#' \item \code{beta_hat} Graph fused lasso estimates
#' \item \code{connected_comp} Estimated connected component
#' \item \code{Naive} the naive p-value using a z-test
#' \item \code{Hyun} the p-value proposed in Hyun et al. (2018)
#' \item \code{hyun_set} the conditioning set of  Hyun et al. (2018) stored as \code{Intervals} class
#' \item \code{CI_result} confidence interval of level 1-\code{alpha_level} if \code{compute_ci=TRUE}
#' }
#' @export
#'
#' @details
#' Consider the generative model \eqn{Y_j = \beta_j + \epsilon_j, \epsilon_j \sim N(0, \sigma^2). j=1,...,n}, where
#' the underlying signal \eqn{\beta} is assumed to be piecewise constant with respect to an underlying
#' graph. The fused lasso estimate minimizes the following objective function
#' \deqn{minimize_{\beta} \frac{1}{2} \sum_{j=1}^{n} ( y_j - \beta_j )^2 + \lambda \sum_{(i,j)\in E}|\beta_i-\beta_j|,}
#' where E is the edge set of the underlying graph. The solution \eqn{\hat{\beta}} can then be
#' segment into connected components; that is, the set of \eqn{\hat{\beta}} that takes on the
#' same value, and are connected in the original graph.
#'
#' Now suppose we want to test whether the means of two estimated connected components \code{c1} and \code{c2}
#' are equal; or equivalently, the null hypothesis of the form \eqn{H_{0}:  \nu^T \beta = 0} versus
#' \eqn{H_{1}:  \nu^T \beta \neq 0} for suitably chosen \eqn{\nu}.
#'
#' This function computes the following p-value:
#' \deqn{P(|\nu^T Y| \ge |\nu^T y| \; | \;  \hat{C}_1, \hat{C}_2 \in CC_K(Y),  \Pi_\nu^\perp Y  =  \Pi_\nu^\perp y),}
#' where \eqn{CC_K(Y)} is the set of estimated connected components from applying K steps of the dual path algorithm on data Y
#' , and \eqn{\Pi_\nu^\perp} is the orthogonal projection to the orthogonal complement of \eqn{\nu}.
#' In particular, the test based on this p-value controls the selective Type I error and has higher power than an existing method
#' by Hyun et al. (2018). Readers can refer to the Section 3 in Chen et al. (2021+) for more details.
#'
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
#' result_demo <- fusedlasso_inf(y=y, D=Dmat, c1=1, c2=2, method="K", sigma=sigma, K=K)
#' summary(result_demo)
#'
#' This tests the null hypothesis of no difference in means between clusters k1 and k2
#' at level K in a hierarchical clustering. (The K clusters are numbered as per the
#' results of the cutree function in the stats package.)
#' @references
#' Chen YT, Witten DM. (2022+) Selective inference for k-means clustering. arXiv preprint.
#' https://arxiv.org/abs/xxxx.xxxxx.
#' @export
kmeans_inference <- function(X, k, cluster_1, cluster_2, sig=NULL, SigInv=NULL,
                             iter.max = 10, seed = 1234, nstart = 1,
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
  if((min(cluster_1,cluster_2)<1)|(max(cluster_1,cluster_2)>k)){
    stop("Cluster numbers must be between 1 and k!")
  }
  n <- dim(X)[1]
  p <- dim(X)[2]
  # get the list of all assigned clusters first
  estimated_k_means <- kmeans_estimation(X, k, iter.max, nstart, seed, tol_eps, verbose)
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



  #final_interval <- intervals::Intervals(c(-Inf,Inf))

  # loop through the initialization
  #init_list <- estimated_k_means$random_init_obs
  #init_cluster <- all_T_clusters[1,]
  # look at covariance matrices -- a different S needed to be computed
  #for (i in c(1:n)){
  #  current_j_prime <- init_cluster[i]
  #  j_star_quad <- norm_sq_phi(XTv, XTv_norm, dir_XTv, v_norm, i,init_list[current_j_prime])
  #  for (j in c(1:length(init_list))){
  #    current_j_quad <- norm_sq_phi(XTv, XTv_norm, dir_XTv, v_norm,i,init_list[j])
  #    curr_quad <- minus_quad_ineq(j_star_quad, current_j_quad)
  #    curr_interval <- intervals::Intervals(solve_one_ineq(curr_quad$quad,
  #                                                         curr_quad$linear, curr_quad$constant))
  #    final_interval <- intervals::interval_intersection(final_interval, curr_interval)

  #  }
  #}

  # loop through all sequence t
  #if(T_length>1){
  #for (l in c(1:(T_length-1))){
   # current_cl <- all_T_clusters[(l+1),]
  #  last_cl <- all_T_clusters[(l),]
    # loop through all the observations
   # for (i in c(1:n)){
      # loop through all cluster classes
    #  current_cl_i <- current_cl[i]
      #
     # k_star_quad <- norm_phi_canonical_kmeans(XTv, XTv_norm, dir_XTv, v_norm, last_cl, current_cl_i, v_vec, i) #i is the observation
    #  for (j in c(1:k)){

     #   k_current_quad <- norm_phi_canonical_kmeans(XTv, XTv_norm, dir_XTv, v_norm, last_cl, j, v_vec, i) #i is the observation
    #    curr_quad <- minus_quad_ineq(k_star_quad, k_current_quad)
    #    curr_interval <- intervals::Intervals(solve_one_ineq(curr_quad$quad,
     #                                                        curr_quad$linear, curr_quad$constant))
        # interval update
     #   final_interval <- intervals::interval_intersection(final_interval, curr_interval)

      #}
  #  }
  #}
  #}

  # final intervals look correct -- try to find truncation?
  #final_interval_chisq <- intervals::interval_intersection(intervals::Intervals(c(0,Inf)),
  #                                                         final_interval)

  # compute test_stats in the isotropic case
  #if(!is.null(sig)){
    #test_stats <- norm_vec(diff_means)
    #scale_factor <- squared_norm_nu*sig^2
  #  gestat <- intervals::Intervals(c(test_stats^2/scale_factor, Inf))
  #}

  #if(!is.null(SigInv)){
    #test_stats <- sqrt(as.numeric(t(diff_means)%*%SigInv%*%diff_means))
    #scale_factor <- squared_norm_nu
   #gestat <- intervals::Intervals(c(test_stats^2/scale_factor, Inf))
  #}

  # improve numerical stability
  final_interval_chisq <- intervals::interval_union(final_interval_chisq,
                           intervals::Intervals_full(c(test_stats-(1e-09),test_stats+(1e-09))),
                           check_valid=FALSE)


  # if (max(as.matrix(final_interval_chisq)[,2])<test_stats){
  #   if(all.equal(max(as.matrix(final_interval_chisq)[,2]),test_stats)){
  #     final_interval_chisq <- intervals::interval_union(final_interval_chisq,
  #                              intervals::Intervals(c(max(final_interval_chisq[,2])
  #                                                     ,test_stats+(1e-09))))
  #   }
  # }

  denom <- final_interval_chisq^2/scale_factor

  #cat("gestat",gestat,"\n")
  #cat("denom",denom,"\n")

  numer <- suppressWarnings(intervals::interval_intersection(gestat, denom))
  pval <- TChisqRatioApprox(p, numer, denom)

  result_list <- list("estimated_k_means" = estimated_k_means, "final_interval"=final_interval_chisq,
                      "final_cluster" = estimated_final_cluster, "test_stats"=test_stats,
                      "sig" = sig, "scale_factor" = scale_factor, "pval" = pval)
  return(result_list)

}







