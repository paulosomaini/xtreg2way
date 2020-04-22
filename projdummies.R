# projdummies R dev
# Brian Chivers, DARC Team

projdummies <- function(hhid, tid, w){
  library("Matrix")
  ## NAN and INF check for w
  # end function if so
  if(any(is.infinite(w))){
    print("Infinite values in w not supported")
    return()
  }
  if(any(is.nan(w))){
    print("NAN values in w not supported")
    return()
  }
  ######################
  hhid_fac <- as.factor(hhid)
  tid_fac <- as.factor(tid)
  
  
  obs <- nrow(w)
  N <- nlevels(hhid_fac)
  T <- nlevels(tid_fac)
  
  DH <- sparseMatrix(i=as.numeric(hhid_fac),j=as.numeric(tid_fac),x=w,dims=list(N,T))
  DD <- rowSums(DH)
  HH <- colSums(DH)
  
  DH <- DH[,1:dim(DH)[2]-1]
  HH <- HH[1:(length(HH)-1)]
  
  invHH <- sparseMatrix(i=1:(T-1), j=1:(T-1), x=HH^-1, dims=list(T-1,T-1))
  invDD <- sparseMatrix(i=1:N, j=1:N, x=DD^-1, dims=list(N,N))
  
    return_list <- list()

  if (N<T){
      A <- solve(diag(DD) - DH %*% invHH %*% t(DH))
      invHHDH <- invHH %*% t(DH)
      B <- -A %*% DH %*% invHH
      
      return_list$A <- A
      return_list$invHHDH <- invHHDH
      return_list$invHH <- invHH
      return_list$B <- B
    } else {
      #Done
      # Complex Conj() missing
      C <- solve(diag(HH)- (t(DH) %*% invDD %*% DH), sparse=TRUE)
      invDDDH=invDD %*% DH
      B= -invDD %*% DH %*% C
      
      return_list$B <- B
      return_list$C <- C
      return_list$invDD <- invDD
      return_list$invDDDH <- invDDDH
    }
    return(return_list)
}
