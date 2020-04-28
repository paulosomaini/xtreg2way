# avar R dev
# Brian Chivers, DARC Team

avar <- function(X, e, group, J){
  library("pracma")
  L <- dim(X)[1]
  K <- dim(X)[2]
  L2 <- dim(e)[1]
  
  ## Check if dimensions match
  if (L != L2){
    print("Infinite values in w not supported")
    #return()
  }
  
  # If J isn't passed, construct it
  if (!exists("J")) {
    J <- t(X) %*% X
  }
  
  G <- nlevels(group)
  
  eX <- sparseMatrix(i=1:L, j=1:L, x=as.numeric(e), dims=list(L,L)) %*% X
  
  if ((length(group)!=L) | (G==L) ){
    V <- t(eX) %*% eX
  } else {
    eX <- sparseMatrix(i=as.numeric(group), j=1:L, x=1, dims=list(G,L)) %*% eX
    V <- t(eX) %*% eX
  }
  
  matCov <- mrdivide(mldivide(J,as.matrix(V)),J)
  return(matCov)
}
