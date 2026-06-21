#' More powerful test for a difference in means between k-means clusters
#'
#' This function computes the original Lloyd-path selective p-value and a
#' union selective p-value that conditions only on preserving the two tested
#' observed clusters along the one-dimensional perturbation path.
#'
#' @param X Numeric matrix; \eqn{n} by \eqn{q} matrix of observed data.
#' @param k Integer; the number of clusters for k-means clustering.
#' @param cluster_1,cluster_2 Two different integers in \eqn{\{1,\ldots,k\}}.
#' @param iso Boolean. If \code{TRUE}, an isotropic covariance model is used with
#' \code{sig}. If \code{FALSE}, supply \code{SigInv} for the general-covariance
#' test.
#' @param sig Numeric; noise standard deviation, used when \code{iso=TRUE}. If
#' omitted, a robust median-based estimator is used, as in
#' \code{\link{kmeans_inference}}.
#' @param SigInv Numeric \eqn{q} by \eqn{q} matrix specifying \eqn{\Sigma^{-1}};
#' required when \code{iso=FALSE}. The Lloyd-path exploration is unchanged (it is
#' Euclidean); only the test statistic and the interval scaling differ, exactly
#' as in \code{\link{kmeans_inference}}.
#' @param iter.max Positive integer; maximum Lloyd iterations.
#' @param seed Random seed for the shared k-means initialization.
#' @param tol_eps A small number specifying the convergence criterion.
#' @param verbose Boolean; currently reserved for compatibility.
#' @param probe_eps Positive number used when probing just outside discovered
#' interval endpoints.
#' @param endpoint_tol Positive number used to decide whether two interval
#' endpoints are numerically identical.
#' @param max_paths Maximum number of distinct Lloyd paths to explore.
#' @param max_halving Maximum number of probe halvings near an endpoint.
#' @param scan_tail_prob Positive tail probability used to bound how far along
#' the perturbation path we keep exploring distinct Lloyd paths. Scanning stops
#' once \eqn{\phi} exceeds the point where the conditional \eqn{\chi} survival
#' mass falls below \code{scan_tail_prob}; regions beyond contribute negligibly
#' to the p-value. Mirrors the truncation used in the generalized-lasso union
#' test and prevents excessive iteration with no practical effect.
#' @param phi_max Optional positive number giving an explicit upper bound on the
#' perturbation parameter \eqn{\phi} to scan. Overrides \code{scan_tail_prob}
#' when supplied.
#'
#' @return A list with the original observed-path interval and p-value
#' (\code{interval_path}, \code{pval_path}), the union interval and p-value
#' (\code{interval_union}, \code{pval_union}), and a long-form \code{regions}
#' table with one row per interval component.
#' @export
kmeans_inference_union <- function(X, k, cluster_1, cluster_2,
                                   iso=TRUE, sig=NULL, SigInv=NULL,
                                   iter.max = 10, seed = 1234,
                                   tol_eps = 1e-6, verbose=TRUE,
                                   probe_eps = 1e-6,
                                   endpoint_tol = 1e-7,
                                   max_paths = 250,
                                   max_halving = 30,
                                   scan_tail_prob = 1e-6,
                                   phi_max = NULL){
  if(!is.matrix(X)) stop("X should be a matrix")
  if(sum(is.na(X))>0){stop("NA is not allowed in the input data X")}
  if(k>=nrow(X)){
    stop("Cannot have more clusters than observations")
  }
  if ((iso)&(is.null(sig))){
    cat("Variance not specified, using a robust median-based estimator by default!\n")
    sig <- .kmeans_estimate_MED(X)
  }
  if(is.null(sig)&is.null(SigInv)){
    stop("At least one of variance and covariance matrix must be specified!")
  }
  if((!is.null(sig))&(!is.null(SigInv))){
    stop("Only one of variance and covariance matrix can be specified!")
  }
  if((!iso)&(is.null(SigInv))){
    stop("You must specify SigInv when iso=FALSE!")
  }
  if(!is.null(SigInv)){
    if(!is.matrix(SigInv) || nrow(SigInv)!=ncol(X) || ncol(SigInv)!=ncol(X)){
      stop("SigInv must be a q by q matrix matching ncol(X).")
    }
  }
  if((min(cluster_1,cluster_2)<1)|(max(cluster_1,cluster_2)>k)){
    stop("Cluster numbers must be between 1 and k!")
  }

  n <- dim(X)[1]
  p <- dim(X)[2]
  estimated_k_means <- kmeans_estimation(X, k, iter.max, seed, tol_eps, verbose)
  if(length(unique(estimated_k_means$final_cluster))<k){
    stop("k-means clustering did not return the desired number of clusters! Try a different seed?")
  }
  estimated_final_cluster <- estimated_k_means$cluster[[estimated_k_means$iter]]

  v_vec <- rep(0, times=nrow(X))
  v_vec[estimated_final_cluster == cluster_1] <- 1/(sum(estimated_final_cluster == cluster_1))
  v_vec[estimated_final_cluster == cluster_2] <- -1/(sum(estimated_final_cluster == cluster_2))

  n1 <- sum(estimated_final_cluster == cluster_1)
  n2 <- sum(estimated_final_cluster == cluster_2)
  squared_norm_nu <- 1/n1 + 1/n2
  v_norm <- sqrt(squared_norm_nu)
  diff_means <- colMeans(X[estimated_final_cluster == cluster_1, ,drop=FALSE]) -
    colMeans(X[estimated_final_cluster == cluster_2, , drop=FALSE])
  XTv <- diff_means
  XTv_norm <- norm_vec(diff_means)
  if(XTv_norm <= 0){
    stop("The observed difference in means has zero norm; the perturbation direction is undefined.")
  }
  dir_XTv <- XTv/XTv_norm

  # The Lloyd-path exploration below is purely Euclidean (it perturbs X along
  # dir_XTv and re-runs Euclidean k-means), so it is identical for the isotropic
  # and general-covariance cases. The covariance only changes the test statistic
  # and a linear reparameterization of the perturbation line: along the fixed 1-D
  # direction, the general-cov statistic sqrt(delta' SigInv delta) equals the
  # Euclidean phi times stat_scale = Sig_XTv_norm / XTv_norm. We therefore
  # explore in the Euclidean phi scale and rescale the resulting intervals by
  # stat_scale before computing the p-value. This matches kmeans_compute_S_genCov,
  # which is kmeans_compute_S_iso with coefficients scaled by 1/stat_scale.
  if(!is.null(sig)){
    test_stat_used    <- XTv_norm
    scale_factor_used <- squared_norm_nu*sig^2
    stat_scale        <- 1
    p_naive <- multivariate_Z_test(X, estimated_final_cluster,
                                   cluster_1, cluster_2, sig)$pval
  } else {
    Sig_XTv_norm <- sqrt(as.numeric(t(diff_means)%*%SigInv%*%diff_means))
    if(Sig_XTv_norm <= 0){
      stop("The observed difference in means has zero SigInv-norm; the test is undefined.")
    }
    test_stat_used    <- Sig_XTv_norm
    scale_factor_used <- squared_norm_nu
    stat_scale        <- Sig_XTv_norm/XTv_norm
    p_naive <- pchisq(test_stat_used^2/scale_factor_used, df=p, lower.tail=FALSE)
  }

  # Bound how far along the perturbation path we scan for distinct Lloyd paths.
  # The chi-ratio p-value integrates a chi density over the conditioning set,
  # with numerator restricted to phi >= test_stat. We scan until the survival
  # mass beyond phi_max is at most scan_tail_prob of the mass beyond the observed
  # statistic, so dropping the far tail changes the p-value by at most that
  # relative amount. This places phi_max strictly above the statistic and scales
  # with the effect size, mirroring the truncation used in the gen-lasso union.
  # phi_max is computed in the statistic scale, then mapped back to the Euclidean
  # exploration scale via stat_scale.
  if(is.null(phi_max)){
    if(!is.numeric(scan_tail_prob) || scan_tail_prob <= 0 || scan_tail_prob >= 1){
      stop("scan_tail_prob must be in (0, 1).")
    }
    stat_chisq   <- test_stat_used^2/scale_factor_used
    surv_stat    <- pchisq(stat_chisq, df=p, lower.tail=FALSE)
    surv_target  <- max(scan_tail_prob*surv_stat, .Machine$double.xmin)
    phi_max_stat <- sqrt(scale_factor_used*qchisq(surv_target, df=p, lower.tail=FALSE))
    if(!is.finite(phi_max_stat)){
      phi_max_stat <- test_stat_used + 10*sqrt(scale_factor_used)  # backstop
    }
    phi_max <- phi_max_stat/stat_scale   # back to Euclidean exploration scale
  }
  if(!is.numeric(phi_max) || length(phi_max) != 1 || phi_max <= 0){
    stop("phi_max must be a single positive number.")
  }
  phi_max <- max(phi_max, XTv_norm)

  path_exploration <- .kmeans_explore_union_regions(
    X=X, k=k, observed_fit=estimated_k_means,
    observed_cluster=estimated_final_cluster,
    cluster_1=cluster_1, cluster_2=cluster_2,
    XTv=XTv, XTv_norm=XTv_norm, dir_XTv=dir_XTv,
    v_vec=v_vec, v_norm=v_norm, iter.max=iter.max,
    tol_eps=tol_eps, probe_eps=probe_eps,
    endpoint_tol=endpoint_tol, max_paths=max_paths,
    max_halving=max_halving, phi_max=phi_max
  )

  # Regions are discovered in the Euclidean perturbation scale; rescale to the
  # statistic scale (identity when isotropic) so the intervals and p-values are
  # expressed on the same scale as the test statistic, matching kmeans_inference.
  interval_path <- .kmeans_scale_interval(.kmeans_regions_to_interval(
    path_exploration$regions[
      path_exploration$regions$preserves_observed_path, , drop=FALSE]), stat_scale)
  interval_union <- .kmeans_scale_interval(.kmeans_regions_to_interval(
    path_exploration$regions[
      path_exploration$regions$preserves_tested_clusters, , drop=FALSE]), stat_scale)

  pval_path <- .kmeans_interval_pvalue(interval_path, test_stat_used, scale_factor_used, p)
  pval_union <- .kmeans_interval_pvalue(interval_union, test_stat_used, scale_factor_used, p)

  result_list <- list("interval_union"=interval_union,
                      "interval_path"=interval_path,
                      "final_interval"=interval_union,
                      "final_cluster" = estimated_final_cluster,
                      "test_stat"=test_stat_used,
                      "cluster_1" = cluster_1,
                      "cluster_2" = cluster_2,
                      "sig" = sig, "SigInv" = SigInv,
                      "scale_factor" = scale_factor_used,
                      "p_naive" = p_naive,
                      "pval_union" = pval_union,
                      "pval_path" = pval_path,
                      "pval" = pval_union,
                      "regions" = path_exploration$regions,
                      "paths" = path_exploration$paths,
                      "exhaustive" = path_exploration$exhaustive,
                      "call" = match.call())
  class(result_list) <- c("kmeans_inference_union", "kmeans_inference")
  return(result_list)
}

.kmeans_estimate_MED <- function(X){
  for (j in c(1:ncol(X))){
    X[,j] <- X[,j]-median(X[,j])
  }
  sqrt(median(X^2)/qchisq(1/2,df=1))
}

.kmeans_fast_dist_compute <- function(x,y) {
  outer(rowSums(x^2), rowSums(y^2), '+') - tcrossprod(x, 2 * y)
}

.kmeans_estimation_fixed_init <- function(X, k, initial_sample, iter.max = 10,
                                          tol_eps = 1e-4, verbose=TRUE){
  if(!is.matrix(X)) stop("X should be a matrix")
  if(k>=nrow(X)){
    stop("Cannot have more clusters than observations")
  }
  iter_T <- 0
  cluster_assign_list <- vector("list", length = iter.max)
  centroid_list <- vector("list", length = iter.max)
  objective_value <- vector("list", length = iter.max)
  current_centroid <- X[initial_sample,,drop=FALSE]
  distance_matrix <- .kmeans_fast_dist_compute(current_centroid, X)
  current_cluster <- apply(distance_matrix,2,which.min)
  if(length(unique(current_cluster)) < k){
    stop("k-means produced an empty cluster along the perturbation path")
  }
  iter_T <- iter_T+1
  centroid_list[[iter_T]] <- current_centroid
  cluster_assign_list[[iter_T]] <- current_cluster
  curr_objective_value <- sum(apply(distance_matrix,2,min))
  objective_value[[iter_T]] <- curr_objective_value
  same_cluster <- FALSE
  while((iter_T<=iter.max)&(!same_cluster)){
    for (current_k in c(1:k)){
      X_current <- X[(current_cluster==current_k), ,drop=FALSE]
      if(nrow(X_current) == 0){
        stop("k-means produced an empty cluster along the perturbation path")
      }
      current_centroid[current_k,] <- .colMeans(X_current, nrow(X_current), ncol(X_current))
    }
    distance_matrix <- .kmeans_fast_dist_compute(current_centroid,X)
    current_cluster <- apply(distance_matrix,2,which.min)
    if(length(unique(current_cluster)) < k){
      stop("k-means produced an empty cluster along the perturbation path")
    }
    iter_T <- iter_T+1
    centroid_list[[iter_T]] <- current_centroid
    cluster_assign_list[[iter_T]] <- current_cluster
    same_cluster <- all(current_cluster==cluster_assign_list[[iter_T-1]])
    new_objective_value <- sum(apply(distance_matrix,2,min))
    curr_objective_value <- new_objective_value
    objective_value[[iter_T]] <- curr_objective_value
  }
  list("cluster" = cluster_assign_list, "centers" = centroid_list,
       "objective" = objective_value, "iter" = iter_T,
       "final_cluster" = cluster_assign_list[[iter_T]],
       "random_init_obs" = initial_sample)
}

.kmeans_perturb_data <- function(X, phi, XTv_norm, dir_XTv, v_vec, v_norm){
  X + as.numeric((phi - XTv_norm)/(v_norm^2)) * tcrossprod(v_vec, dir_XTv)
}

.kmeans_path_signature <- function(fit){
  cluster_list <- fit$cluster[seq_len(fit$iter)]
  paste(vapply(cluster_list, function(z) paste(z, collapse=","), character(1)),
        collapse="|")
}

.kmeans_centroids_from_path <- function(X, cluster_list, initial_sample, k){
  centers <- vector("list", length(cluster_list))
  centers[[1]] <- X[initial_sample,,drop=FALSE]
  if(length(cluster_list) > 1){
    for(idx in 2:length(cluster_list)){
      last_cl <- cluster_list[[idx-1]]
      curr_centers <- matrix(NA_real_, nrow=k, ncol=ncol(X))
      for(current_k in c(1:k)){
        X_current <- X[(last_cl==current_k), ,drop=FALSE]
        if(nrow(X_current) == 0){
          stop("Cannot compute a path region with an empty previous cluster")
        }
        curr_centers[current_k,] <- .colMeans(X_current, nrow(X_current), ncol(X_current))
      }
      centers[[idx]] <- curr_centers
    }
  }
  centers
}

.kmeans_compute_path_interval_iso <- function(X, fit, XTv, XTv_norm, dir_XTv,
                                             v_vec, v_norm, k){
  cluster_list <- fit$cluster[seq_len(fit$iter)]
  all_T_clusters <- do.call(rbind, cluster_list)
  all_T_centroids <- .kmeans_centroids_from_path(
    X, cluster_list, fit$random_init_obs, k
  )
  interval <- kmeans_compute_S_iso(X, fit, all_T_clusters, all_T_centroids,
                                   nrow(X), XTv, XTv_norm, dir_XTv,
                                   v_vec, v_norm, nrow(all_T_clusters), k)
  .kmeans_interval_positive(interval)
}

.kmeans_interval_empty <- function(){
  intervals::Intervals_full(matrix(numeric(0), ncol=2))
}

.kmeans_interval_positive <- function(interval){
  if(nrow(interval) == 0) return(.kmeans_interval_empty())
  positive <- suppressWarnings(intervals::interval_intersection(
    interval, intervals::Intervals_full(c(0, Inf))
  ))
  if(nrow(positive) == 0) return(.kmeans_interval_empty())
  intervals::reduce(positive, check_valid=FALSE)
}

.kmeans_scale_interval <- function(interval, s){
  if(s == 1 || nrow(interval) == 0) return(interval)
  if(!is.finite(s) || s <= 0) stop("interval scale factor must be finite and positive")
  intervals::reduce(intervals::Intervals_full(as.matrix(interval)*s),
                    check_valid=FALSE)
}

.kmeans_interval_from_matrix <- function(mat){
  if(is.null(mat) || nrow(mat) == 0) return(.kmeans_interval_empty())
  intervals::reduce(intervals::Intervals_full(mat), check_valid=FALSE)
}

.kmeans_interval_contains <- function(interval, phi, tol=1e-8){
  if(!is.finite(phi) || nrow(interval) == 0) return(FALSE)
  mat <- as.matrix(interval)
  any(phi >= mat[,1] - tol & phi <= mat[,2] + tol)
}

.kmeans_interval_component_containing <- function(interval, phi, tol=1e-8){
  if(!is.finite(phi) || nrow(interval) == 0) return(NULL)
  mat <- as.matrix(interval)
  idx <- which(phi >= mat[,1] - tol & phi <= mat[,2] + tol)
  if(length(idx) == 0) return(NULL)
  mat[idx[1],]
}

.kmeans_interval_to_region_rows <- function(interval, path_id, anchor_phi,
                                           preserves_observed_path,
                                           preserves_tested_clusters){
  mat <- as.matrix(interval)
  if(nrow(mat) == 0){
    return(data.frame(region_id=integer(), set_id=character(),
                      path_id=character(), lower=numeric(), upper=numeric(),
                      anchor_phi=numeric(),
                      preserves_observed_path=logical(),
                      preserves_tested_clusters=logical(),
                      stringsAsFactors=FALSE))
  }
  data.frame(region_id=seq_len(nrow(mat)), set_id=path_id,
             path_id=path_id, lower=mat[,1], upper=mat[,2],
             anchor_phi=anchor_phi,
             preserves_observed_path=preserves_observed_path,
             preserves_tested_clusters=preserves_tested_clusters,
             stringsAsFactors=FALSE)
}

.kmeans_regions_to_interval <- function(regions){
  if(is.null(regions) || nrow(regions) == 0){
    return(.kmeans_interval_empty())
  }
  .kmeans_interval_from_matrix(as.matrix(regions[,c("lower", "upper"), drop=FALSE]))
}

.kmeans_add_region_rows <- function(regions, new_rows){
  if(nrow(new_rows) == 0) return(regions)
  start <- nrow(regions)
  new_rows$region_id <- start + seq_len(nrow(new_rows))
  rbind(regions, new_rows)
}

.kmeans_preserves_observed_clusters <- function(observed_cluster, candidate_cluster,
                                                cluster_1, cluster_2){
  preserves_one <- function(cluster_id){
    obs_idx <- which(observed_cluster == cluster_id)
    candidate_labels <- unique(candidate_cluster[obs_idx])
    if(length(candidate_labels) != 1) return(FALSE)
    all(observed_cluster[candidate_cluster == candidate_labels] == cluster_id)
  }
  preserves_one(cluster_1) && preserves_one(cluster_2)
}

.kmeans_add_endpoint_probes <- function(probe_queue, probe_seen, interval,
                                       endpoint_tol, phi_max=Inf){
  if(nrow(interval) == 0) return(probe_queue)
  mat <- as.matrix(interval)
  for(i in seq_len(nrow(mat))){
    lower <- mat[i,1]
    upper <- mat[i,2]
    if(is.finite(lower) && lower > endpoint_tol && lower <= phi_max){
      key <- paste0("L:", signif(lower, 14))
      if(is.null(probe_seen[[key]])){
        probe_seen[[key]] <- TRUE
        probe_queue <- rbind(probe_queue,
                             data.frame(endpoint=lower, direction=-1))
      }
    }
    # Only probe rightward past an endpoint that is still below phi_max; beyond
    # phi_max any further Lloyd paths contribute negligible chi mass.
    if(is.finite(upper) && upper < phi_max){
      key <- paste0("R:", signif(upper, 14))
      if(is.null(probe_seen[[key]])){
        probe_seen[[key]] <- TRUE
        probe_queue <- rbind(probe_queue,
                             data.frame(endpoint=upper, direction=1))
      }
    }
  }
  probe_queue
}

.kmeans_explore_union_regions <- function(X, k, observed_fit, observed_cluster,
                                          cluster_1, cluster_2,
                                          XTv, XTv_norm, dir_XTv, v_vec, v_norm,
                                          iter.max, tol_eps, probe_eps,
                                          endpoint_tol, max_paths, max_halving,
                                          phi_max=Inf){
  observed_signature <- .kmeans_path_signature(observed_fit)
  initial_sample <- observed_fit$random_init_obs
  paths <- list()
  signature_to_path <- new.env(parent=emptyenv())
  probe_seen <- new.env(parent=emptyenv())
  regions <- data.frame(region_id=integer(), set_id=character(),
                        path_id=character(), lower=numeric(), upper=numeric(),
                        anchor_phi=numeric(),
                        preserves_observed_path=logical(),
                        preserves_tested_clusters=logical(),
                        stringsAsFactors=FALSE)
  probe_queue <- data.frame(endpoint=numeric(), direction=integer())
  exhaustive <- TRUE

  register_anchor <- function(anchor_phi){
    if(!is.finite(anchor_phi) || anchor_phi < 0 || anchor_phi > phi_max) return(NULL)
    discovered <- .kmeans_regions_to_interval(regions)
    if(.kmeans_interval_contains(discovered, anchor_phi, endpoint_tol)){
      return(NULL)
    }
    X_anchor <- .kmeans_perturb_data(X, anchor_phi, XTv_norm, dir_XTv, v_vec, v_norm)
    fit <- tryCatch(
      .kmeans_estimation_fixed_init(X_anchor, k, initial_sample, iter.max,
                                    tol_eps=tol_eps, verbose=FALSE),
      error=function(e) NULL
    )
    if(is.null(fit)) return(NULL)
    signature <- .kmeans_path_signature(fit)
    if(!is.null(signature_to_path[[signature]])){
      return(signature_to_path[[signature]])
    }

    path_id <- paste0("path_", sprintf("%03d", length(paths)+1))
    signature_to_path[[signature]] <- path_id
    interval <- .kmeans_compute_path_interval_iso(
      X, fit, XTv, XTv_norm, dir_XTv, v_vec, v_norm, k
    )
    preserves_observed_path <- identical(signature, observed_signature)
    preserves_tested_clusters <- .kmeans_preserves_observed_clusters(
      observed_cluster, fit$final_cluster, cluster_1, cluster_2
    )
    paths[[path_id]] <<- list(path_id=path_id,
                              signature=signature,
                              anchor_phi=anchor_phi,
                              interval=interval,
                              cluster_path=fit$cluster[seq_len(fit$iter)],
                              final_cluster=fit$final_cluster,
                              iter=fit$iter,
                              preserves_observed_path=preserves_observed_path,
                              preserves_tested_clusters=preserves_tested_clusters)
    rows <- .kmeans_interval_to_region_rows(
      interval, path_id, anchor_phi,
      preserves_observed_path, preserves_tested_clusters
    )
    regions <<- .kmeans_add_region_rows(regions, rows)
    probe_queue <<- .kmeans_add_endpoint_probes(
      probe_queue, probe_seen, interval, endpoint_tol, phi_max
    )
    path_id
  }

  register_anchor(XTv_norm)

  while(nrow(probe_queue) > 0){
    if(length(paths) >= max_paths){
      exhaustive <- FALSE
      warning("Reached max_paths before exhausting the perturbation path.")
      break
    }
    probe <- probe_queue[1,,drop=FALSE]
    probe_queue <- probe_queue[-1,,drop=FALSE]
    endpoint <- probe$endpoint
    direction <- probe$direction
    step <- max(probe_eps, abs(endpoint)*probe_eps)
    accepted <- FALSE
    for(halve in 0:max_halving){
      candidate <- endpoint + direction*step
      if(direction < 0 && candidate < 0){
        candidate <- endpoint/2
      }
      if(candidate < 0 || !is.finite(candidate)){
        break
      }
      # Reaching the scan boundary is a deliberate stop, not a failure to
      # verify adjacency: regions beyond phi_max are intentionally skipped.
      if(direction > 0 && candidate > phi_max){
        accepted <- TRUE
        break
      }
      discovered <- .kmeans_regions_to_interval(regions)
      if(.kmeans_interval_contains(discovered, candidate, endpoint_tol)){
        accepted <- TRUE
        break
      }
      path_id <- register_anchor(candidate)
      if(is.null(path_id)){
        step <- step/2
        next
      }
      component <- .kmeans_interval_component_containing(
        paths[[path_id]]$interval, candidate, endpoint_tol
      )
      if(is.null(component)){
        step <- step/2
        next
      }
      adjoining_endpoint <- if(direction > 0) component[1] else component[2]
      if(abs(adjoining_endpoint - endpoint) <= endpoint_tol + abs(endpoint)*endpoint_tol){
        accepted <- TRUE
        break
      }
      step <- step/2
    }
    if(!accepted){
      exhaustive <- FALSE
      warning("Could not verify adjacency for at least one perturbation endpoint.")
    }
  }

  list(regions=regions, paths=paths, exhaustive=exhaustive)
}

.kmeans_interval_pvalue <- function(interval, test_stats, scale_factor, df){
  if(nrow(interval) == 0) return(NA_real_)
  interval <- intervals::interval_union(
    interval,
    intervals::Intervals_full(c(test_stats-(1e-09), test_stats+(1e-09))),
    check_valid=FALSE
  )
  denom <- interval^2/scale_factor
  gestat <- intervals::Intervals(c(test_stats^2/scale_factor, Inf))
  numer <- suppressWarnings(intervals::interval_intersection(gestat, denom))
  if(nrow(numer) == 0) return(0)
  TChisqRatioApprox(df, numer, denom)
}

#' @export
summary.kmeans_inference_union <- function(object, ...){
  data.frame(cluster_1 = object$cluster_1,
             cluster_2 = object$cluster_2,
             test_stat = object$test_stat,
             p_union = object$pval_union,
             p_path = object$pval_path,
             p_naive = object$p_naive,
             exhaustive = object$exhaustive)
}
