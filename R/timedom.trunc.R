#' This function removes lags from a linear filter.
#'
#' @title Choose lags of an object of class \code{timedom}
#' @param A an object of class \code{timedom}.
#' @param lags a vector which contains a set of lags. These lags must be a subset of the lags defined for timedom object A. Only those lags will be kept, the other lags are removed.
#' @return An object of class \code{timedom}.
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

