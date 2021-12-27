# ----- functions for solving quadratic inequalities -----
# ----- adopted from Lucy Gao's function here:  -----
# -----  https://github.com/lucylgao/clusterpval/blob/master/R/trunc_sets.R -----
#' Solve the roots of quadratic polynomials related to testing for a difference in means
#'
#' Solves \eqn{ax^2 + bx + c \leq 0}. If the solution is empty, the function returns NA.
#'
#' @keywords internal
#'
#' @param A, B, C the coefficients of the quadratic equation.
#' @param tol if \eqn{|a|}, \eqn{|b|}, or \eqn{|c|} is not larger than tol, then treat it as zero.
#'
#' @return Returns an "Intervals" object containing the solution set.
solve_one_ineq <- function(A, B, C, tol=1e-8) {
  # Computes the complement of the set {phi: B*phi + C <=  0},
  compute_linear_ineq_complement <- function(B, C, tol=1e-8) {
    #  If B = 0
    if(abs(B) <= tol) {
      if(C <= tol) { # C <= 0: inequality is always satisfied
        return(c(-Inf,Inf)) # all of real line
      } else { # C > 0: something has gone wrong -- no solution works
        warning("B = 0 and C > 0: B*phi + C <= 0 has no solution")
        return(c(0,0)) # do not return any value
      }
    }

    # If B \neq 0
    ratio <- -C/B
    # If B > 0:
    if(B > tol) {
      return(c(-Inf,ratio))
      }
    if(B < tol) {
      return(c(ratio, Inf))
    }
    # we will return a degenerate interval
    return(c(0,0))
  }


  # A = 0?
  if(abs(A) <= tol) {
    return(compute_linear_ineq_complement(B, C, tol))
  }

  # We know A \neq 0
  discrim <- B^2 - 4*A*C

  # If discriminant is small, we assume there is no root
  if(discrim <= tol) {
    if(A > tol) { # Parabola opens up: there is no solution
      return(c(0,0))
    } else { # Parabola opens down: every x is a solution
      return(c(-Inf, Inf))
    }
  }

  # We now know that A =/= 0, and that there are two roots
  # we compute the roots using the suggestion outlined at
  # https://people.csail.mit.edu/bkph/articles/Quadratics.pdf
  # for numerical stability purposes
  sqrt_discrim <- sqrt(discrim)
  if (B >= tol){
    root_1 <- (-B-sqrt_discrim)/(2*A)
    root_2 <- (2*C)/(-B-sqrt_discrim)
    roots <- sort(c(root_1, root_2))
  }else{
    root_1 <- (-B+sqrt_discrim)/(2*A)
    root_2 <- (2*C)/(-B+sqrt_discrim)
    roots <- sort(c(root_1, root_2))
  }

  # Parabola opens up? (A > 0?)
  if(A > tol) {
    return(roots)
  }else{
    interval_result <- matrix(c(-Inf,roots[1], roots[2], Inf),
           ncol=2, byrow = T)
    return(interval_result)
  }

  # if everything fails -- then we do not have a solution
  warning("Edge case for quadratic inequality solver!")
  return(c(0,0))
}

# ----- functions for solving quadratic inequalities -----
#' Solve the roots of quadratic polynomials related to testing for a difference in means
#'
#' Solves \eqn{ax^2 + bx + c \ge 0}, then returns the complement of the solution set
#' wrt to the real line, unless the complement is empty, in which case
#' the function returns NA.
#'
#' @keywords internal
#'
#' @param A, B, C the coefficients of the quadratic equation.
#' @param tol if \eqn{|a|}, \eqn{|b|}, or \eqn{|c|} is not larger than tol, then treat it as zero.
#'
#' @return Returns an "Intervals" object containing NA or the complement of the solution set.
solve_one_ineq_complement <- function(A, B, C, tol=1e-10) {
  # Computes the complement of the set {phi: B*phi + C <=  0},
  compute_linear_ineq_complement <- function(B, C, tol=1e-8) {
    #  If B = 0
    if(abs(B) <= tol) {
      if(C <= tol) { # C <= 0: inequality is always satisfied
        return(c(0,0)) # all of real line
      } else { # C > 0: something has gone wrong -- no solution works
        warning("B = 0 and C > 0: B*phi + C <= 0 has no solution")
        return(c(-Inf,Inf)) # do not return any value
      }
    }

    # If B \neq 0
    ratio <- -C/B
    # If B > 0:
    if(B > tol) {
      return(c(ratio,Inf))
    }
    if(B < tol) {
      return(c(-Inf, ratio))
    }
    # we will return a degenerate interval
    # return(c(-Inf,Inf))
  }


  # A = 0?
  if(abs(A) <= tol) {
    return(compute_linear_ineq_complement(B, C, tol))
  }

  # We know A \neq 0
  discrim <- B^2 - 4*A*C

  # If discriminant is small, we assume there is no root
  if(discrim <= tol) {
    if(A > tol) { # Parabola opens up: there is no solution
      return(c(-Inf,Inf))
    } else { # Parabola opens down: every x is a solution
      return(c(0, 0))
    }
  }

  # We now know that A =/= 0, and that there are two roots
  # we compute the roots using the suggestion outlined at
  # https://people.csail.mit.edu/bkph/articles/Quadratics.pdf
  # for numerical stability purposes
  sqrt_discrim <- sqrt(discrim)
  if (B >= tol){
    root_1 <- (-B-sqrt_discrim)/(2*A)
    root_2 <- (2*C)/(-B-sqrt_discrim)
    roots <- sort(c(root_1, root_2))
  }else{
    root_1 <- (-B+sqrt_discrim)/(2*A)
    root_2 <- (2*C)/(-B+sqrt_discrim)
    roots <- sort(c(root_1, root_2))
  }

  if(A > tol) {
    if(roots[1] > tol) {
      return(c(0, roots[1], roots[2], Inf))
    }

    if(roots[2] <= tol) {
      warning("something wrong with the discriminant calculation!")
      return(c(0,Inf))
    }

    return(c(roots[2], Inf))
  }

  # We now know that there are two roots, and parabola opens down (A < 0)
  if(roots[2] < -tol) {
    return(c(0, 0))
  }

  if(roots[1] > tol) {
    return(c(roots[1], roots[2]))
  }

  return(c(-Inf, roots[2]))

  # Parabola opens up? (A > 0?)
  # if(A > tol) {
  #
  #   if(roots[1] > tol){
  #     return()
  #   }
  #
  #   interval_result <- matrix(c(-Inf,roots[1], roots[2], Inf),
  #                             ncol=2, byrow = T)
  #
  #   return(interval_result)
  # }else{
  #   return(roots)
  # }

  # if everything fails -- then we do not have a solution
  #warning("Edge case for quadratic inequality solver!")
  #return(c(-Inf,Inf))
}


#' Represent <x'(phi)_i, x'(phi)_j> as a quadratic function in phi ----
#' @keywords internal
#'
#' @param X, matrix n by p
#' @param v, contrast vector n by 1
#' @param i, first index
#' @param j, second index
#'
#' @return parameters: a, b, c the coefficients of the quadratic equation such that (ax^2 + bx + c <= 0)
#'
inner_product_phi <- function(X, v, i, j){
  v_norm <- norm_vec(v)
  XTv <- t(X)%*%v
  XTv_norm <- norm_vec(XTv)
  dir_XTv <- XTv/norm_vec(XTv)
  quad_coef <- (v[i]*v[j])/(v_norm)^4
  linear_coef <- (v[j]/v_norm^2)*(X[i,]%*%dir_XTv) + (v[i]/v_norm^2)*(X[j,]%*%dir_XTv) -
    2*v[i]*v[j]*XTv_norm/(v_norm^4)
  constant_vec_1 <- X[i,]-v[i]*XTv_norm/v_norm^2*dir_XTv
  constant_vec_2 <- X[j,]-v[j]*XTv_norm/v_norm^2*dir_XTv
  constant_coef <- sum(constant_vec_1*constant_vec_2)
  coef_list <- list("quad" = as.numeric(quad_coef), "linear" = as.numeric(linear_coef),
                    "constant"= as.numeric(constant_coef))
  return(coef_list)
}


#' Represent ||x'(phi)_i-x'(phi)_j||_2^2 as a quadratic function in phi ----
#' @keywords internal
#'
#' @param XTv, vector p by 1
#' @param XTv_norm, norm of XTv
#' @param dir_XTv, vector p by 1 := XTv/XTv_norm
#' @param v_norm, 2-norm of vector v
#' @param i, first index
#' @param j, second index
#' @param v, the vector
#'
#' @return parameters: a, b, c the coefficients of the quadratic equation
#' such that (ax^2 + bx + c <= 0)
#'
norm_sq_phi <- function(X, v, XTv, XTv_norm, dir_XTv, v_norm, i, j){
  #v_norm <- norm_vec(v)
  #XTv <- t(X)%*%v
  #XTv_norm <- norm_vec(XTv)
  #dir_XTv <- XTv/norm_vec(XTv)
  quad_coef <- (v[i]-v[j])^2/(v_norm)^4
  linear_coef <- 2*((((v[i]-v[j])/v_norm^2)*(X[i,]-X[j,])%*%dir_XTv) - ((v[i]-v[j])/v_norm^2)^2*XTv_norm)
  constant_vec <- X[i,]-X[j,]-(v[i]-v[j])*XTv/(v_norm^2)
  constant_coef <- sum(constant_vec*constant_vec)
  coef_list <- list("quad" = as.numeric(quad_coef), "linear" = as.numeric(linear_coef),
                    "constant"= as.numeric(constant_coef))
  return(coef_list)
}




#' Represent <x'(phi)_i, x'(phi)_j> as a quadratic function in phi ----
#' @keywords internal
#'
#' @param XTv, vector p by 1
#' @param XTv_norm, norm of XTv
#' @param dir_XTv, vector p by 1 := XTv/XTv_norm
#' @param v_norm, 2-norm of vector v
#' @param cl, factor vector n by 1 (most recent cluster assignment)
#' @param k, cluster of interest
#' @param v, contrast vector n by 1
#' @param i, index of observation
#'
#' @return parameters: a, b, c the coefficients of the quadratic equation such that (ax^2 + bx + c <= 0)
#'
norm_phi_canonical_kmeans <- function(X, last_centroids, XTv, XTv_norm, dir_XTv, v_norm, cl, k, v, i){
  n_k <- sum(cl==k)
  indicator_vec <- rep(0, times=length(cl))
  indicator_vec[cl==k] <- 1
  indicator_location <- which(cl==k)
  #v_norm <- norm_vec(v)
  #XTv <- t(X)%*%v
  #XTv_norm <- norm_vec(XTv)
  #dir_XTv <- XTv/norm_vec(XTv)
  # compute quad coef
  v_i_expression <- (v[i]-sum(indicator_vec*v)/n_k)/(v_norm^2)
  #X_current <- X[indicator_location,]
  #x_i_expression <- X[i,] - ((t(indicator_vec) %*% X)/n_k)
  x_i_expression <- X[i,] - last_centroids[k,]
  # n_k and class k
  # .colMeans(X[indicator_location,], n_k, dim(X)[2])
  quad_coef <- (v_i_expression)^2
  # compute lienar coef
  linear_coef_part_1 <- v_i_expression*(x_i_expression%*%dir_XTv)
  linear_coef_part_2 <- (v_i_expression)^2*XTv_norm
  linear_coef <- 2*(linear_coef_part_1-linear_coef_part_2)

  constant_vec <- x_i_expression - c(v_i_expression*XTv)
  constant_coef <- sum(constant_vec*constant_vec)

  coef_list <- list("quad" = as.numeric(quad_coef), "linear" = as.numeric(linear_coef),
                    "constant"= as.numeric(constant_coef))
  return(coef_list)
}

#' Implement the minus operation for two quadratic inequalities
#'
minus_quad_ineq <- function(quad1, quad2){
  coef_list <- list("quad" = quad1$quad-quad2$quad, "linear" = quad1$linear-quad2$linear,
                    "constant"= quad1$constant-quad2$constant)
  return(coef_list)
}


#' Compute set S for isotropic case
#' @export
kmeans_compute_S_iso <- function(X, estimated_k_means, all_T_clusters,
                                 all_T_centroids,
                                 n, XTv, XTv_norm,
                                 dir_XTv, v_vec,v_norm,T_length, k){

  final_interval <- intervals::Intervals(c(0,Inf))
  all_interval_lists <- list()

  # loop through the initialization
  init_list <- estimated_k_means$random_init_obs
  init_cluster <- all_T_clusters[1,]
  # look at covariance matrices -- a different S needed to be computed
  for (i in c(1:n)){
    current_j_prime <- init_cluster[i]
    j_star_quad <- norm_sq_phi(X, v_vec, XTv, XTv_norm, dir_XTv, v_norm,i,init_list[current_j_prime])
    for (j in c(1:length(init_list))){
      current_j_quad <- norm_sq_phi(X, v_vec, XTv, XTv_norm, dir_XTv, v_norm,i,init_list[j])
      curr_quad <- minus_quad_ineq(j_star_quad, current_j_quad)
      curr_interval <- solve_one_ineq_complement(curr_quad$quad, curr_quad$linear, curr_quad$constant)
      all_interval_lists[[(i-1)*length(init_list)+j]] <- curr_interval
      #final_interval <- intervals::interval_intersection(final_interval, curr_interval)

    }
  }

  curr_len <- length(all_interval_lists)
  curr_counter <- 1
  # loop through all sequence t
  if(T_length>1){
    for (l in c(1:(T_length-1))){
      current_cl <- all_T_clusters[(l+1),]
      last_cl <- all_T_clusters[(l),]
      last_centroids <- all_T_centroids[[(l+1)]] # k by q matrix
      # pre-compute the centroids (or maybe extract from kmeans estimation??)
      # loop through all the observations
      for (i in c(1:n)){
        # loop through all cluster classes
        current_cl_i <- current_cl[i]
        k_star_quad <- norm_phi_canonical_kmeans(X, last_centroids, XTv, XTv_norm, dir_XTv, v_norm, last_cl,
                                                 current_cl_i, v_vec, i) #i is the observation
        for (j in c(1:k)){
          if(j!=current_cl_i){
            k_current_quad <- norm_phi_canonical_kmeans(X, last_centroids, XTv, XTv_norm,
                                                        dir_XTv, v_norm, last_cl, j, v_vec, i) #i is the observation
            curr_quad <- minus_quad_ineq(k_star_quad, k_current_quad)
            curr_interval <- solve_one_ineq_complement(curr_quad$quad,
                                                       curr_quad$linear,
                                                       curr_quad$constant)
            all_interval_lists[[curr_len+curr_counter]] <- curr_interval
            curr_counter <- curr_counter + 1
          }
          # interval update
          #final_interval <- intervals::interval_intersection(final_interval, curr_interval)

        }
      }
    }
  }

  # final intervals look correct -- try to find truncation?
  final_interval_complement <- do.call('c', all_interval_lists)
  final_interval_complement <- matrix(final_interval_complement, ncol=2, byrow=T)
  final_interval_complement <- intervals::reduce(intervals::Intervals_full(final_interval_complement),
                                                     check_valid=FALSE)

  final_interval_chisq <- intervals::interval_complement(final_interval_complement)
  #intervals::interval_intersection(intervals::Intervals(c(0,Inf)),
                          #                                 final_interval)
  return(final_interval_chisq)
}



#' Compute set S for isotropic case
#' @export
kmeans_compute_S_genCov <- function(X, estimated_k_means, all_T_clusters,
                                    all_T_centroids,
                                 n, XTv, XTv_norm,
                                 dir_XTv, v_vec,v_norm,T_length,
                                 Sig_XTv_norm, k){

  Sig_Inv_factor <- XTv_norm/Sig_XTv_norm
  final_interval <- intervals::Intervals(c(0,Inf))
  # keep track of all the intervals
  all_interval_lists <- list()

  # loop through the initialization
  init_list <- estimated_k_means$random_init_obs
  init_cluster <- all_T_clusters[1,]
  # look at covariance matrices -- a different S needed to be computed
  for (i in c(1:n)){
    current_j_prime <- init_cluster[i]
    j_star_quad <- norm_sq_phi(X, v_vec, XTv, XTv_norm, dir_XTv, v_norm,i,init_list[current_j_prime])
    for (j in c(1:length(init_list))){
      current_j_quad <- norm_sq_phi(X, v_vec, XTv, XTv_norm, dir_XTv, v_norm,i,init_list[j])
      curr_quad <- minus_quad_ineq(j_star_quad, current_j_quad)
      curr_interval <- solve_one_ineq_complement((Sig_Inv_factor)^2*curr_quad$quad,
                                                           (Sig_Inv_factor)*curr_quad$linear,
                                                           curr_quad$constant)
      all_interval_lists[[(i-1)*length(init_list)+j]] <- curr_interval
      #final_interval <- intervals::interval_intersection(final_interval, curr_interval)
    }
  }

  # keep track of the list elements
  curr_len <- length(all_interval_lists)
  curr_counter <- 1
  # loop through all sequence t
  if(T_length>1){
    for (l in c(1:(T_length-1))){
      current_cl <- all_T_clusters[(l+1),]
      last_cl <- all_T_clusters[(l),]
      # get pre-computed centroids
      last_centroids <- all_T_centroids[[(l+1)]] # k by q matrix
      # loop through all the observations
      for (i in c(1:n)){
        # loop through all cluster classes
        current_cl_i <- current_cl[i]
        k_star_quad <- norm_phi_canonical_kmeans(X, last_centroids, XTv, XTv_norm, dir_XTv, v_norm, last_cl,
                                                 current_cl_i, v_vec, i) #i is the observation
        for (j in c(1:k)){

          k_current_quad <- norm_phi_canonical_kmeans(X, last_centroids, XTv, XTv_norm,
                                                      dir_XTv, v_norm, last_cl, j, v_vec, i) #i is the observation
          curr_quad <- minus_quad_ineq(k_star_quad, k_current_quad)
          curr_interval <- solve_one_ineq_complement((Sig_Inv_factor)^2*curr_quad$quad,
                                                               (Sig_Inv_factor)*curr_quad$linear,
                                                               curr_quad$constant)
          # interval update
          all_interval_lists[[curr_len+curr_counter]] <- curr_interval
          curr_counter <- curr_counter + 1

        }
      }
    }
  }

  # final intervals look correct -- try to find truncation?
  final_interval_complement <- do.call('c', all_interval_lists)
  final_interval_complement <- matrix(final_interval_complement, ncol=2, byrow=T)
  final_interval_complement <- intervals::reduce(intervals::Intervals(final_interval_complement),
                                                 check_valid=FALSE)

  final_interval_chisq <- intervals::interval_complement(final_interval_complement)
  #intervals::interval_intersection(intervals::Intervals(c(0,Inf)),
  #                                 final_interval)
  return(final_interval_chisq)

}



