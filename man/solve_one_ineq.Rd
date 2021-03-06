% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trunc_sets.R
\name{solve_one_ineq}
\alias{solve_one_ineq}
\title{Solve the roots of quadratic polynomials related to testing for a difference in means}
\usage{
solve_one_ineq(A, B, C, tol = 1e-08)
}
\arguments{
\item{A, }{B, C the coefficients of the quadratic equation.}

\item{tol}{if \eqn{|a|}, \eqn{|b|}, or \eqn{|c|} is not larger than tol, then treat it as zero.}
}
\value{
Returns an "Intervals" object containing the solution set.
}
\description{
Solves \eqn{ax^2 + bx + c \leq 0}. If the solution is empty, the function returns NA.
}
\keyword{internal}
