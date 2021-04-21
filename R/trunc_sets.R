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


