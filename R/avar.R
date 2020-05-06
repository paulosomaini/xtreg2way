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
#' hhid <- c("a","b","c")
#' tid <- c("1","2","3")
#' w <- c(1,1,1)
#' x1 <- c(1,2,3)
#' y <- c(2,4,6)
#' struc <- projdummies(hhid, tid, w)
#' x1_projected <- projvar(x1, struc)
#' y_projected <- projvar(y, struc)
#' reg <- regress1(y_projected, x1_projected)
#' matCov <- avar(X, reg$res, struc$hhid, reg$XX)

avar <- function(X, e, group=NULL, J=NULL) {
  L <- dim(X)[1]
  L2 <- dim(e)[1]

  ## Check if dimensions match
  if (L != L2) {
    print("The 1st dimensions of X and e must match")
    return()
  }

  # If J isn't passed, construct it
  if (is.null(J)) {
    J <- t(as.matrix(X)) %*% X
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

  matCov <- pracma::mrdivide(pracma::mldivide(J, as.matrix(V)), J)
  return(matCov)
}
