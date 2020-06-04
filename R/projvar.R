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
#' @return A list will be returned with the following named values:
#'         var - the projected variable
#'         delta, tau - intermediate variables
#'
#'
#' @examples
#' hhid <- c("a","b","c","a","b","c" ,"a","b","c" ,"a","b","c" ,"a","b","c")
#' tid <- c("1","1" ,"1" ,"2","2" ,"3","3","3" ,"4","4","5" ,"5","6","6" ,"6")
#' w <- rep(1, 15)
#' x1 <- rnorm(15, mean=50, sd=10)   
#' 
#' struc <- projdummies(hhid, tid, w)
#' x1p <- projvar(x1, struc)
#' @export

projvar <- function(var, struc) {
    ## NAN and INF check for w
    # end function if so
    if (any(is.infinite(var))) {
      stop("Infinite values in w not supported")
    }
    if (any(is.nan(var))) {
      stop("NAN values in w not supported")
    }
  
    if ("esample" %in% names(struc) ) {
      if (length(var) == length(struc$esample)) {
        var <- var[struc$esample]
      }
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
    
    return_list <- list()
    return_list$delta <- delta
    return_list$tau <- tau
    return_list$var <- (var - delta[c(struc$hhid)] - tau[c(struc$tid)]) * sqrt(struc$w)

    return(return_list)
}
