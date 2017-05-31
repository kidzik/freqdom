#' This function removes lags from a linear filter.
#'
#' @title Truncates a time-domain object to given lags
#' @param A an object of class \code{timedom}.
#' @param lags a vector which contains a set of lags. These lags must be a subset of the lags defined for timedom object A. Only those lags will be kept, the other lags are removed.
#' @return An object of class \code{timedom}.
#' @examples
#' X = rar(100)
#' Y = rar(100)
#' #estimate regressors in model $Y_t = \sum_{i\in Z} A_i X_{t-i}$
#' A = speclagreg(X,Y)
#' B = timedom.trunc(A, c(-1, 2, 3))
#' @export 
timedom.trunc = function(A, lags){
  if (!is.timedom(A))
    stop ("A must be a time domain filter")
  
  D = dim(A$operators)
  D[1] = sum(A$lags %in% lags)
  A$operators = array(A$operators[A$lags %in% lags,,],D)
  A$lags = intersect(lags,A$lags)
  A
}

