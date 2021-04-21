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






