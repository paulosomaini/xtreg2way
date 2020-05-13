# regress1 R dev
# Brian Chivers, DARC Team

#' Regression
#'
#' \code{regress1} performs an OLS regression based on the projected variables y and X.
#'
#'
#' @param X A matrix or vector of independent variable(s)
#' @param y The dependent variable
#' @return A list which contains X'X, the returned coefficients beta, and residuals res
#'
#'
#' @examples
#' hhid <- c("a","b","c","a","b","c" ,"a","b","c" ,"a","b","c" ,"a","b","c")
#' tid <- c("1","1" ,"1" ,"2","2" ,"3","3","3" ,"4","4","5" ,"5","6","6" ,"6")
#' w <- rep(1, 15)
#' x1 <- rnorm(15, mean=50, sd=10)   
#' x2 <- rnorm(15, mean=50, sd=10)
#' y <- x1 + rnorm(15, mean=50, sd=10)
#' 
#' struc <- projdummies(hhid, tid, w)
#' x1p <- projvar(x1, struc)
#' x2p <- projvar(x2, struc)
#' yp <- projvar(y, struc)
#' 
#' reg <- regress1(yp, data.frame(x1p,x2p))
#' @export

regress1 <- function(y, X) {
  X <- Matrix::Matrix(as.matrix(X))
  XX <- Matrix::t(X) %*% X
  beta <- pracma::mldivide(as.matrix(XX), as.matrix(Matrix::t(X) %*% y))
  res <- y - X %*% beta

  return_list <- list()
  return_list$XX <- XX
  return_list$beta <- beta
  return_list$res <- res
  return(return_list)
}
