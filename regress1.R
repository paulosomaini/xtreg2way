# regress1 R dev
# Brian Chivers, DARC Team
regress1 <- function(y, X){
  library("pracma")
  XX <- t(X) %*% X
  beta <- mldivide(XX, t(X)%*%y)
  res <- y - X %*% beta
  
  return_list <- list()
  return_list$XX <- XX
  return_list$beta <- beta
  return_list$res <- res
  return(return_list)
}

