# regress1 R dev
# Brian Chivers, DARC Team

#' Regression
#'
#' \code{regress1} performs an OLS regression based on the projected variables y and X.
#'
#'
#' @param X A matrix or vector of independant variable(s)
#' @param y The dependant variable
#' @return A list which contains X'X, the returned coefficients beta, and residuals res
#'
#'
#' @examples
#' hhid <- c("a","b","c")
#' tid <- c("1","2","3")
#' w <- c(1,1,1)
#' x1 <- c(1,2,3)
#' y <- c(2,4,6)
#' struc <- projdummies(hhid, tid, w)
#' x1_projected <- projvar(x1, struc)
#' y_projected <- projvar(y, struc)
#' reg <- regress1(y_projected, x1_projected)

regress1 <- function(y, X) {
  XX <- t(X) %*% X
  beta <- pracma::mldivide(XX, t(X) %*% y)
  res <- y - X %*% beta

  return_list <- list()
  return_list$XX <- XX
  return_list$beta <- beta
  return_list$res <- res
  return(return_list)
}
