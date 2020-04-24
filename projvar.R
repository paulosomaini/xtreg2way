#projvar
#Brian Chivers, DARC team

projvar <- function(var, struc){
    ## NAN and INF check for w
    # end function if so
    library("Matrix")
    if(any(is.infinite(var))){
      print("Infinite values in w not supported")
      return()
    }
    if(any(is.nan(var))){
      print("NAN values in w not supported")
      return()
    }
    N <- nlevels(struc$hhid)
    T <- nlevels(struc$tid)
    if (N<T){
      A <- struc$A
      B <- struc$B
      invHH <- Matrix(struc$invHH, sparse=T)
      invHHDH <- struc$invHHDH
    } else{
      B <- struc$B
      C <- struc$C
      invDD <- Matrix(struc$invDD, sparse=T)
      invDDDH <- struc$invDDDH
    }
    
    
    
    aux <- sparseMatrix(i=as.numeric(struc$hhid),j=as.numeric(struc$tid),x=var*struc$w,dims=list(N,T))

    Dy <- rowSums(aux)
    Ty <- colSums(aux)
    Ty <- Ty[1:length(Ty)-1]

    if (N<T){
      delta <- A %*% Dy + B %*% Ty;
      tau <- crossprod(B, Dy) + invHH %*% Ty + invHHDH %*% A %*% crossprod(invHHDH,Ty)
      tau <- rbind(tau,0)
    } else{
      delta <- invDD %*% Dy + invDDDH %*% C %*% crossprod(invDDDH,Dy) + B %*% Ty
      tau <- crossprod(B,Dy) + C %*% Ty
      tau <- rbind(tau,0)

    }

    delta <- as.numeric(delta)
    tau <- as.numeric(tau)

    return((var-delta[c(struc$hhid)]-tau[c(struc$tid)])*sqrt(struc$w))
}




