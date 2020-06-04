# avar R dev
# Brian Chivers, DARC Team

#' Asymptotic variance of Estimator
#'
#' \code{avar} calculated the asymptotic variance of the regression estimation
#'
#'
#' @param X A matrix or vector of independent variable(s)
#' @param e The residuals from the regression
#' @param group (optional) The cluster identifier (hhid from \code{projdummies})
#' @param J (optional) This is assumed to be X'X, and can be input if pre-calculated
#' @return A matrix of the covariates
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
#' projvar_list <- projvar(x1, struc)
#' x1p <- projvar_list$var
#' projvar_list <- projvar(x2, struc)
#' x2p <- projvar_list$var
#' projvar_list <- projvar(y, struc)
#' yp <- projvar_list$var
#' 
#' reg <- regress1(yp, data.frame(x1p,x2p))
#' matCov <- avar(data.frame(x1p, x2p), reg$res, struc$hhid, reg$XX)
#' @export

avar <- function(X, e, group=NULL, J=NULL) {
  X <- Matrix::Matrix(as.matrix(X))
  L <- dim(X)[1]
  L2 <- dim(e)[1]

  ## Check if dimensions match
  if (L != L2) {
    print("The 1st dimensions of X and e must match")
    return()
  }

  # If J isn't passed, construct it
  if (is.null(J)) {
    J <- Matrix::t(Matrix::Matrix(X)) %*% X
  }

  G <- nlevels(group)
  eX <- Matrix::sparseMatrix(i = 1:L, j = 1:L,
                     x = as.numeric(e), dims = list(L, L)) %*% X

  if (is.null(group) | (length(group) != L) | (G == L)) {
    V <- Matrix::t(eX) %*% eX
  } else {
    eX <- Matrix::sparseMatrix(i = as.numeric(group), j = 1:L,
                       x = 1, dims = list(G, L)) %*% eX
    V <- Matrix::t(eX) %*% eX
  }

  matCov <- pracma::mrdivide(pracma::mldivide(as.matrix(J), as.matrix(V)), as.matrix(J))
  return(matCov)
}
