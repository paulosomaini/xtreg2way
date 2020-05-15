# xtreg2way
# Brian Chivers, DARC Team

#' Estimation of Two Way Fixed Effects Model
#'
#' \code{xtreg2way} estimates a 2-way fixed effect model absorbing the two set of dummies and reports standard errors.
#'
#'
#' @param y The dependant variable, size N-by-1
#' @param X The matrix of covariates, size N-by-K
#' @param iid The group ID, size N-by-1
#' @param tid The time ID, size N-by-1
#' @param w (optional) The vector of weights, size N-by-1.  If omitted, w will be set to 1 for all observations
#' @param struc (optional) This list contains the results from the first step of the algorithm.  To save computational time, you can rerun the algorithm on different columns by providing this struc.
#' @param se (optional) This indicates standard error estimate to be calulcated.  Possible values include:
#'        se=="0" : standard errors assuming homoscedasticity and no within  group correlation or serial correlation.
#'        se=="1" : standard errors  proposed by Arellano (1987) robust to heteroscedasticity and serial correlation.
#'        se=="2" : standard errors robust to heteroscedasticity but assumes no correlation within group or serial correlation.
#'        se=="11" : Arellano standard errors with a degree of freedom correction performed by Stata xtreg, fe.
#'        If se is omitted or set to [] then it is set to 1 and the Arellano (1987) estimator is computed.
#' @param noise (optional) If noise is set to "1", then results are displayed
#' @param data (optional) If the first argument is a formula, you can pass a dataframe here
#' @return betaHat (K-by-1) vector of estimated coefficients
#' aVarHat (K-by-K) estimate of the matrix of variances and covariance of  the estimator.
#' yp (N-by-1) the residual of the projection of y on the two sets of  dummies.
#' Xp (N-by-K) the residual of the projection of each column of X on the two  sets of dummies.
#' struc (structure) results of the first step of the algorithm.

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
#' output3 <- xtreg2way(y~x1+x2, iid=hhid, tid=tid, w=w, se="2", noise="1",
#'                     data=data.frame(x1,x2,y))
#' @export


xtreg2way <- function(y, X=NULL, iid = NULL, tid = NULL, w = NULL, 
                      struc = NULL, se = "", noise = NULL, data = NULL) {
  #This variable is needed for the return at the bottom
  struc_is_null <- is.null(struc)
  #If struc is passed, grab iid tid and w from it
  if (!is.null(struc)) {
    iid <- struc$hhid
    tid <- struc$tid
    w <- struc$w
  }
  #If a formula is passed, grab all necessary variables
  if( class(y) == "formula") {
    
    #These variables should be passed explicitly or via struc by now
    if (is.null(iid) | is.null(tid)) {
      stop(paste("Arguments Missing:",
                 "If the first argument is a formula,",
                 "then iid and tid are required arguments"), sep = " ")
    }
    
    #Check to see if labels exist
    for (label in attr(stats::terms(y),"term.labels")) {
      if (!label %in% labels(data)[[2]]){
        stop(paste("Missing Variable:",
                   label,"is not in data", sep=" "))
      }
    }
    #Build X
    X <- data[attr(stats::terms(y),"term.labels")]
    
    #Make sure only 1 y variable exists
    y_label <- setdiff(all.vars(y) , attr(stats::terms(y),"term.labels"))
    if( length(y_label) > 1) {
      stop(paste("Multiple independent variables in formula",
                 y_label, sep=" "))
    }
    y <- data[y_label]
  }
  
  X <- as.matrix(X)
  y <- as.matrix(y)
  obs <- dim(X)[1]
  K <- dim(X)[2]

  #If w is null, fill with ones
  if (is.null(w)) {
    w <- rep(1, obs)
  }

  #checking y and X dimensions
  if (length(y) != obs) {
    stop(paste("Error, Dimension Mismatch:",
          "y must be a vector of length N,",
          "and X must be a N by K matrix"), sep = " ")
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

    #Project variables
    for (col_id in 1:K) {
      X[, col_id] <- projvar(X[, col_id], struc)
    }
    y <- projvar(y, struc)

  } else {
    #Else, if struc is provided, check length of hhid and tid
    if (length(struc$hhid) != obs | length(struc$tid) != obs) {
      stop(paste("Error, Dimension Mismatch:",
            "struc$hhid and struc$tid must be vectors of length N,",
            "matching the length of y and rows of X"), sep = " ")
    }
  }

  #Perform regression on projected variables
  reg <- regress1(y, X)
  betaHat <- Matrix::t(reg$beta)

  #SE == '0' for standard errors
  #assuming homoscedasticity and no within group correlation
  #or serial correlation
  if (se == "0") {
    N <- nlevels(struc$hhid)
    T <- nlevels(struc$tid)
    sig2hat <- (Matrix::t(reg$res) %*% reg$res) /
      (sum(struc$w > 0) - N - T + 1 - length(reg$beta))
    aVarHat <- sqrt(diag(as.numeric(sig2hat) * Matrix::solve(reg$XX)))
  #SE=='1'
  #standard errors proposed by Arellano (1987) robust to
  #heteroscedasticity and serial correlation
  } else if (se == "1") {
    aVarHat <- avar(X, reg$res, struc$hhid, reg$XX)
  #SE == 2
  #it computes standard errors robust to heteroscedasticity,
  #but assumes no correlation within group or serial correlation.
  } else if (se == "2") {
    aVarHat <- avar(X, reg$res, as.factor(1:obs), reg$XX)
  #SE == 11
  #Arellano (1987) standard errors with a degree of freedom
  #correction performed by Stata xtreg, fe
  } else if (se == "11") {
    aVarHat <- avar(X, reg$res, struc$hhid, reg$XX)
    N <- nlevels(struc$hhid)
    stata_dof <- ((obs - 1) / (obs - length(reg$beta) - 1)) * (N / (N - 1))
    aVarHat <- aVarHat * (stata_dof)^2;
  #ELSE
  #Arellano (1987) estimator is computed
  } else {
    aVarHat <- avar(X, reg$res, struc$hhid, reg$XX)
  }

  #Build the noise DF and print it
  if (!is.null(noise)) {
    std <- sqrt(diag(aVarHat))
    df <- data.frame(coefficients = as.numeric(betaHat), se = std,
                     tstat = Matrix::t(betaHat) / std,
                     pval = (1 - stats::pnorm(abs(t(betaHat) / std), 0, 1)) / 2)
    print(df)
  }

  #Return all that is needed
  return_list <- list()
  return_list$betaHat <- betaHat
  return_list$aVarHat <- aVarHat

  if (struc_is_null) {
    return_list$y <- y
    return_list$X <- X
    return_list$struc <- struc
  }
  return(return_list)
}
