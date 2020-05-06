#projvar
#Brian Chivers, DARC team

#' Project Variables
#'
#' \code{projvar} uses the matrices from \code{projdummies} to project a variable.
#' In practice, all dependent and indepdent variables must be project for the two way regression
#'
#'
#' @param var A vector of a single variable
#' @param struc The output of \code{projdummies}, containing matrices necessary to project
#' @return A vector the same length as var will be returned, with each item corresponding to a projected value from var
#'
#'
#' @examples
#' hhid <- c("a","b","c")
#' tid <- c("1","2","3")
#' w <- c(1,1,1)
#' x1 <- c(1,2,3)
#' struc <- projdummies(hhid, tid, w)
#' x1_projected <- projvar(x1, struc)

projvar <- function(var, struc) {
    ## NAN and INF check for w
    # end function if so
    if (any(is.infinite(var))) {
      stop("Infinite values in w not supported")
    }
    if (any(is.nan(var))) {
      stop("NAN values in w not supported")
    }
    N <- nlevels(struc$hhid)
    T <- nlevels(struc$tid)
    if (N < T) {
      A <- struc$A
      B <- struc$B
      invHH <- Matrix::Matrix(struc$invHH, sparse = T)
      invHHDH <- struc$invHHDH
    } else{
      B <- struc$B
      C <- struc$C
      invDD <- Matrix::Matrix(struc$invDD, sparse = T)
      invDDDH <- struc$invDDDH
    }

    aux <- Matrix::sparseMatrix(i = as.numeric(struc$hhid), j = as.numeric(struc$tid),
                        x = as.numeric(var * struc$w), dims = list(N, T))

    Dy <- Matrix::rowSums(aux)
    Ty <- Matrix::colSums(aux)
    Ty <- Ty[1:(length(Ty)-1)]
    if (N < T) {

      delta <- A %*% Dy + B %*% Ty;
      tau <- Matrix::crossprod(B, Dy) +
             invHH %*% Ty +
             invHHDH %*% A %*% Matrix::crossprod(invHHDH, Ty)
      tau <- rbind(tau, 0)
    } else {
      delta <- invDD %*% Dy +
               invDDDH %*% C %*% Matrix::crossprod(invDDDH, Dy) +
               B %*% Ty
      tau <- Matrix::crossprod(B, Dy) + C %*% Ty
      tau <- rbind(tau, 0)

    }

    delta <- as.numeric(delta)
    tau <- as.numeric(tau)

    return((var - delta[c(struc$hhid)] - tau[c(struc$tid)]) * sqrt(struc$w))
}
