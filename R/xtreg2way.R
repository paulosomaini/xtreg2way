# xtreg2way
# Brian Chivers, DARC Team

#' Estimation of Two Way Fixed Effects Model
#'
#' \code{xtreg2way} estimates a 2-way fixed effect model absorbing the two set of dummies and reports standard errors.
#' @aliases xtreg2way.formula
#' @aliases xtreg2way.default
#' 
#' 
#' @param y Either a numeric of data, or a formula
#' @param formula Model specifications
#' @param data A dataframe with labels from the formula \code{y}
#' @param X A matrix of data
#' @param iid (optional) The group ID, size N-by-1 - not needed if \code{struc} is provided
#' @param tid (optional) The time ID, size N-by-1 - not needed if \code{struc} is provided
#' @param w (optional) The vector of weights, size N-by-1.  If omitted, w will be set to 1 for all observations
#' @param struc (optional) This list contains the results from the first step of the algorithm.  To save computational time, you can rerun the algorithm on different columns by providing this struc.
#' @param se (optional) This indicates standard error estimate to be calculcated.  Possible values include:
#'        se=="0" : standard errors assuming homoscedasticity and no within  group correlation or serial correlation.
#'        se=="1" : standard errors  proposed by Arellano (1987) robust to heteroscedasticity and serial correlation.
#'        se=="2" : standard errors robust to heteroscedasticity but assumes no correlation within group or serial correlation.
#'        se=="11" : Arellano standard errors with a degree of freedom correction performed by Stata xtreg, fe.
#'        If se is omitted or set to [] then it is set to 1 and the Arellano (1987) estimator is computed.
#' @param noise (optional) If noise is set to "1", then results are displayed
#' @param ... Other parameters, based on method used
#' @return \code{betaHat} (K-by-1) vector of estimated coefficients
#' 
#' \code{aVarHat} (K-by-K) estimate of the matrix of variances and covariance of  the estimator.
#' 
#' \code{y} (N-by-1) the residual of the projection of y on the two sets of  dummies.
#' 
#' \code{X} (N-by-K) the residual of the projection of each column of X on the two  sets of dummies.
#' 
#' \code{struc} (list) results of the first step of the algorithm.

#'
#'
#' @examples
#' hhid <- c("a","b","c","a","b","c" ,"a","b","c" ,"a","b","c" ,"a","b","c")
#' tid <- c("1","1" ,"1" ,"2","2" ,"3","3","3" ,"4","4","5" ,"5","6","6" ,"6")
#' w <- rep(1, 15)
#' x1 <- rnorm(15, mean=50, sd=10)   
#' x2 <- rnorm(15, mean=50, sd=10)   
#' y <- x1 + rnorm(15, mean=50, sd=10)    
#' #The most basic way to use this function
#' output <- xtreg2way(y, x1, hhid, tid, w, se="2", noise="1")
#' #You can rerun faster with different columns using output$struc
#' output2 <- xtreg2way(y, data.frame(x1,x2), struc=output$struc)
#' #Or you can use a formula and specify data=
#' output3 <- xtreg2way(y~x1+x2, data=data.frame(x1,x2,y), iid=hhid, tid=tid, w=w, 
#'                      se="2", noise="1")
#'            
#' @export

xtreg2way <- function(y, ...){
  UseMethod("xtreg2way",y)
}

#' @describeIn xtreg2way This function ingests a formula as the first argument,
#'  and requires \code{data} as a data.frame
#' @export
xtreg2way.formula <- function(formula, data, iid = NULL, tid = NULL, w = NULL, 
                              struc = NULL, se = "", cluster=NULL, noise = "", ...) {
  # This function inputs a formula
  # and creates variables compatable with xtreg2way.default
  #Check to see if labels exist
  for (label in attr(stats::terms(formula),"term.labels")) {
    if (!label %in% labels(data)[[2]]){
      stop(paste("Missing Variable:",
                 label,"is not in data", sep=" "))
    }
  }
  #Build X
  X <- data[attr(stats::terms(formula),"term.labels")]
  
  #Make sure only 1 y variable exists
  y_label <- setdiff(all.vars(formula) , attr(stats::terms(formula),"term.labels"))
  if( length(y_label) > 1) {
    stop(paste("Multiple independent variables in formula",
               y_label, sep=" "))
  }
  y <- as.matrix(data[y_label])
  
  if(is.character(iid) & length(iid) == 1) {
    if(iid %in% colnames(data)) {
      iid <- data[iid]
    } else {
      stop("If iid is a string, it needs to be a column in data")
    }
  }
  if(is.character(tid) & length(tid) == 1) {
    if(tid %in% colnames(data)){
      tid <- data[tid]
    } else {
      stop("If tid is a string, it needs to be a column in data")
    }
  }
  if(is.character(w) & length(w) == 1) {
    if(w %in% colnames(data)){
      w <- data[w]
    } else {
      stop("If w is a string, it needs to be a column in data")
    }
  }
  
  #NA Checks for X,y
  #The checks for iid, tid, and w are in xtreg2way.default
  if(any(is.na(y))) {
    stop(paste("Error: NA values in the y argument (",
               y_label,")",sep=" "))
  }
  if(any(is.na(X))) {
    stop(paste("Error: NA values in the X argument (",
               attr(stats::terms(formula),"term.labels"),")",sep=" "))
  }

  xtreg2way.default(y, X, iid, tid, w, struc, se, cluster, noise)
}

#' @describeIn xtreg2way Default Method
#' @export
xtreg2way.default<- function(y, X, iid = NULL, tid = NULL, w = NULL, 
                             struc = NULL, se = "", cluster=NULL, noise = "",...) {
  #This variable is needed for the return at the bottom
  struc_is_null <- is.null(struc)
  #If struc is passed, grab iid tid and w from it
  if (!is.null(struc)) {
    iid <- struc$hhid
    tid <- struc$tid
    w <- struc$w
  }
  
  X <- as.matrix(X)
  obs <- dim(X)[1]
  K <- dim(X)[2]
  
  #If w is null, fill with ones
  if (is.null(w)) {
    w <- rep(1, obs)
  }
  
  #NA Checks for X,y, iid, tid, and w
  if(any(is.na(y))) {
    stop("Error: NA values in the y argument")
  }
  
  if(any(is.na(X))) {
    stop("Error: NA values in the X argument")
  }
  if(any(is.na(iid))) {
    stop("Error: NA values in the iid argument")
  }
  if(any(is.na(tid))) {
    stop("Error: NA values in the tid argument")
  }
  if(any(is.na(w))) {
    stop("Error: NA values in the w argument")
  }
  
  #checking y and X dimensions
  if (length(y) != obs) {
    stop(paste("Error, Dimension Mismatch:",
               "y must be a vector of length N,",
               "and X must be a N by K matrix"), sep = " ")
  }
  
  
  ##Redundant Check
  redundant <- xtreg2way::nonredundant(iid, tid, w)
  nr <- redundant$nr
  if (redundant$flag) {
    esample <- (iid %in% nr$iid) & (tid %in% nr$tid)
    y <- y[esample]
    X <- X[esample,]
    iid <- iid[esample]
    tid <- tid[esample]
    w <- w[esample]

    X <- as.matrix(X)
    obs <- dim(X)[1]
    K <- dim(X)[2]
  }
  
  #If struc isn't provided
  if (is.null(struc)) {
    #we need iid and tid
    if (is.null(tid) | is.null(iid)) {
      stop(paste("Error: if struc isn't provided,",
                 "iid and tid are required arguments"))
    }
    #Check length of iid and tid
    if (length(iid) != obs | length(tid) != obs) {
      stop(paste("Error, Dimension Mismatch:",
                 "iid and tid must be vectors of length N,",
                 "matching the length of y and rows of X"), sep = " ")
    }
    #Build struc
    struc <- projdummies(iid, tid, w)
    
    if (redundant$flag) {
      struc$esample <- esample
    }
    
  } else {
    #Else, if struc is provided, check length of hhid and tid
    if (length(struc$hhid) != obs | length(struc$tid) != obs) {
      stop(paste("Error, Dimension Mismatch:",
                 "struc$hhid and struc$tid must be vectors of length N,",
                 "matching the length of y and rows of X"), sep = " ")
    }
  }
  
  #Project variables
  for (col_id in 1:K) {
    projvar_list <- projvar(X[, col_id], struc)
    X[, col_id] <- projvar_list$var
  }
  projvar_list <- projvar(y, struc)
  y <- projvar_list$var
  
  #Perform regression on projected variables
  reg <- regress1(y, X)
  betaHat <- Matrix::t(reg$beta)
  dof <- obs / (obs - length(unique(iid)) - length(unique(tid)) - length(reg$beta)+struc$correction_rank)

  #cluster option
  if (is.null(cluster)){
    cluster<-struc$hhid
  }
  else{
    cluster<-as.factor(cluster)
  }
  
  #SE == '0' for standard errors
  #assuming homoscedasticity and no within group correlation
  #or serial correlation
  if (se == "0") {
    N <- nlevels(struc$hhid)
    T <- nlevels(struc$tid)
    sig2hat <- (Matrix::t(reg$res) %*% reg$res) /
    (sum(struc$w > 0) - N - T + 1 - length(reg$beta)+struc$correction_rank)
    aVarHat <- (as.numeric(sig2hat) * MASS::ginv(as.matrix(reg$XX)))
      (sum(struc$w > 0) - N - T + 1 - length(reg$beta))
    aVarHat <- (as.numeric(sig2hat) * MASS::ginv((as.matrix(reg$XX))))
    #SE=='1'
    #standard errors proposed by Arellano (1987) robust to
    #heteroscedasticity and serial correlation
  } else if (se == "1") {
    aVarHat <- avar(X, reg$res, cluster, reg$XX) * dof
    #SE == 2
    #it computes standard errors robust to heteroscedasticity,
    #but assumes no correlation within group or serial correlation.
  } else if (se == "2") {
    aVarHat <- avar(X, reg$res, as.factor(1:obs), reg$XX) * dof
    #SE == 11
    #Arellano (1987) standard errors with a degree of freedom
    #correction performed by Stata xtreg, fe
  } else if (se == "11") {
    aVarHat <- avar(X, reg$res, cluster, reg$XX)
    N <- nlevels(struc$hhid)
    stata_dof <- ((obs - 1) / (obs - length(reg$beta) - 1)) * (N / (N - 1))
    aVarHat <- aVarHat * (stata_dof)^2;
    #ELSE
    #Arellano (1987) estimator is computed
  } else {
      aVarHat <- avar(X, reg$res, struc$hhid, reg$XX)* dof
  }
  

  #Build the noise DF and print it
  if (noise=="1") {
    std <- sqrt(diag(aVarHat))
    df <- data.frame(coefficients = as.numeric(betaHat), se = std,
                     tstat = Matrix::t(betaHat) / std,
                     pval = (1 - stats::pnorm(abs(t(betaHat) / std), 0, 1)) / 2)
    colnames(df) <- c("coefficients","se","tstat","pval")
    print(df)
  }
  
  
  
  #Return all that is neededz
  return_list <- list()
  return_list$betaHat <- betaHat
  return_list$aVarHat <- aVarHat
  
  if (struc_is_null) {
    return_list$y <- y
    return_list$X <- X
    return_list$struc <- struc
  }
  class(return_list) <- "xtreg2way"
  return(return_list)
}
