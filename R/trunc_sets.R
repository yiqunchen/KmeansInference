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
  # Computes the complement of the set {phi >= 0: B*phi + C >= 0},
  # ignoring (-Inf, 0].
  compute_linear_ineq_complement <- function(B, C, tol=1e-10) {
    # Is B = 0?
    if(abs(B) <= tol) {
      if(C >= -tol) { # C >= 0: inequality automatically satisfied
        return()
      } else { # C < 0: something has gone wrong ...
        warning("B = 0 and C < 0: B*phi + C >=0 is degenerate")
        return(c(0, Inf))
      }
    }

    # We know that B =/= 0
    ratio <- -C/B
    # Is B > 0?
    if(B > tol) {
      if(C >= -tol) { # -C/B <= 0: inequality automatically satisfied
        return()
      } else { # -C/B > 0: the interval extends to the right
        return(c(0, ratio))
      }
    }

    # We know B < 0
    if(C <= tol) { # -C/B <= 0: inequality can't be satisfied
      return(c(0, Inf))
    }

    # We know B < 0 & -C/B > 0: the interval extends to the left
    return(c(ratio, Inf))
  }


  # A = 0?
  if(abs(A) <= tol) {
    return(compute_linear_ineq_complement(B, C, tol))
  }

  # We know A =/= 0
  discrim <- B^2 - 4*A*C

  # No roots or one root?
  if(discrim <= tol) {
    if(A > tol) { # Parabola opens up: inequality automatically satisfied
      return()
    } else { # Parabola opens down: inequality never satisfied
      return(c(0, Inf))
    }
  }

  # We now know that A =/= 0, and that there are two roots
  sqrt_discrim <- sqrt(discrim)
  roots <- sort(c(-B + sqrt_discrim, -B - sqrt_discrim)/(2*A))
  # Parabola opens up? (A > 0?)
  if(A > tol) {
    if(roots[1] > tol) {
      return(c(roots[1], roots[2]))
    }

    if(roots[2] <= tol) {
      return()
    }

    return(c(0, roots[2]))
  }

  # We now know that there are two roots, and parabola opens down (A < 0)
  if(roots[2] < -tol) {
    return(c(0, Inf))
  }

  if(roots[1] > tol) {
    return(c(0, roots[1], roots[2], Inf))
  }

  return(c(roots[2], Inf))
}


### Represent <x'(phi)_i, x'(phi)_j> as a quadratic function in phi ----
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


### Represent ||x'(phi)_i-x'(phi)_j||_2^2 as a quadratic function in phi ----
#' @keywords internal
#'
#' @param X, matrix n by p
#' @param v, contrast vector n by 1
#' @param i, first index
#' @param j, second index
#'
#' @return parameters: a, b, c the coefficients of the quadratic equation such that (ax^2 + bx + c <= 0)
#'
norm_sq_phi <- function(X, v, i,j){
  v_norm <- norm_vec(v)
  XTv <- t(X)%*%v
  XTv_norm <- norm_vec(XTv)
  dir_XTv <- XTv/norm_vec(XTv)
  quad_coef <- (v[i]-v[j])^2/(v_norm)^4
  linear_coef <- 2*((((v[i]-v[j])/v_norm^2)*(X[i,]-X[j,])%*%dir_XTv) - ((v[i]-v[j])/v_norm^2)^2*XTv_norm)
  constant_vec <- X[i,]-X[j,]-(v[i]-v[j])*XTv/(v_norm^2)
  constant_coef <- sum(constant_vec*constant_vec)
  coef_list <- list("quad" = as.numeric(quad_coef), "linear" = as.numeric(linear_coef),
                    "constant"= as.numeric(constant_coef))
  return(coef_list)
}




### Represent <x'(phi)_i, x'(phi)_j> as a quadratic function in phi ----
#' @keywords internal
#'
#' @param X, matrix n by p
#' @param cl, factor vector n by 1 (most recent cluster assignment)
#' @param k, cluster of interest
#' @param v, contrast vector n by 1
#' @param i, index of observation
#'
#' @return parameters: a, b, c the coefficients of the quadratic equation such that (ax^2 + bx + c <= 0)
#'
norm_phi_canonical_kmeans <- function(X, cl, k, v, i){
  n_k <- sum(cl==k)
  indicator_vec <- rep(0, times=length(cl))
  indicator_vec[cl==k] <- 1
  indicator_location <- which(cl==k)
  v_norm <- norm_vec(v)
  XTv <- t(X)%*%v
  XTv_norm <- norm_vec(XTv)
  dir_XTv <- XTv/norm_vec(XTv)
  # compute quad coef
  v_i_expression <- (v[i]-sum(indicator_vec*v)/n_k)/(v_norm^2)
  x_i_expression <- X[i,] - colMeans(X[indicator_location,])
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


combine_quad_ineq <- function(quad1, quad2){
  coef_list <- list("quad" = quad1$quad-quad2$quad, "linear" = quad1$linear-quad2$linear,
                    "constant"= quad1$constant-quad2$constant)
  return(coef_list)
}
