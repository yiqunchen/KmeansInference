#' k means estimation
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


#' k means inference
#' @export
kmeans_inference <- function(X, k, cluster_1, cluster_2, sig=NULL, SigInv=NULL,
                             iter.max = 10, nstart = 1,
                             seed = 1234, tol_eps = 1e-6, verbose=TRUE){

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







