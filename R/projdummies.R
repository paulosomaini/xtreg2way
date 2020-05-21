# projdummies R dev
# Brian Chivers, DARC Team

#' Projection Dummies
#'
#' \code{projdummies} returns necessary matrices to project variables on fixed effect dummies.
#' The input parameters all need to be of the same length.
#'
#'
#' @param hhid A vector of individual effect identifiers
#' @param tid A vector of time effect identifiers
#' @param w A vector of weights for each observation
#' @return A list will be returned with necessary matrices to project upon.
#'    If the time effect has more levels, the matrices B, C, invDD, and invDDDH will be returned
#'    If the individual effect has more levels, the matrices A, B, invHH and invHHDH will be returned
#'
#'    hhid and tid as factors will always be returned, as well as the original weights w that are passed.
#'
#' @examples
#' hhid <- c("a","b","c","a","b","c" ,"a","b","c" ,"a","b","c" ,"a","b","c")
#' tid <- c("1","1" ,"1" ,"2","2" ,"3","3","3" ,"4","4","5" ,"5","6","6" ,"6")
#' w <- rep(1, 15)
#' projdummies(hhid, tid, w)
#' @export


projdummies <- function(hhid, tid, w) {
  ## NAN and INF check for w
  # end function if so
  if (any(is.infinite(w))) {
    stop("Infinite values in w not supported")
  }
  if (any(is.nan(w))) {
    stop("NAN values in w not supported")
  }

  if (length(hhid) != length(tid) | length(hhid) != length(w)) {
    stop("Lengths of hhid, tid, and w must be equal")
  }
  ######################
  hhid_fac <- as.factor(hhid)
  tid_fac <- as.factor(tid)

  return_list <- list()
  return_list$hhid <- hhid_fac
  return_list$tid <- tid_fac
  return_list$w <- w

  N <- nlevels(hhid_fac)
  T <- nlevels(tid_fac)

  DH <- Matrix::sparseMatrix(i = as.numeric(hhid_fac), j = as.numeric(tid_fac),
                     x = w, dims = list(N, T))
  DD <- Matrix::rowSums(DH)
  HH <- Matrix::colSums(DH)

  DH <- DH[, 1:dim(DH)[2] - 1]
  HH <- HH[1:(length(HH) - 1)]

  invHH <- Matrix::sparseMatrix(i = 1 :(T - 1), j = 1 : (T - 1),
                        x = HH^-1, dims = list(T - 1, T - 1))
  invDD <- Matrix::sparseMatrix(i = 1:N, j = 1:N,
                        x = DD^-1, dims = list(N, N))


  if (N < T) {
      A <- MASS::ginv(as.matrix(diag(DD) - DH %*% invHH %*% Matrix::t(DH)))
      invHHDH <- invHH %*% Matrix::t(DH)
      B <- -A %*% DH %*% invHH

      return_list$A <- A
      return_list$invHHDH <- invHHDH
      return_list$invHH <- invHH
      return_list$B <- B
    } else {
      C <- MASS::ginv(as.matrix(diag(HH) - (Matrix::t(DH) %*% invDD %*% DH)))
      invDDDH <- invDD %*% DH
      B <-  -invDD %*% DH %*% C

      return_list$B <- B
      return_list$C <- C
      return_list$invDD <- invDD
      return_list$invDDDH <- invDDDH
    }
    return(return_list)
}
