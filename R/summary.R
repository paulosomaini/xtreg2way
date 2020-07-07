#' @export
summary.xtreg2way <- function (object, ...) {
  std <- sqrt(diag(object$aVarHat))
  df <- data.frame(coefficients = as.numeric(object$betaHat), se = std,
                   tstat = Matrix::t(object$betaHat) / std,
                   pval = (1 - stats::pnorm(abs(t(object$betaHat) / std), 0, 1)) / 2)
  colnames(df) <- c("coefficients","se","tstat","pval")
  print(df)
}