#' This function removes lags from a linear filter.
#'
#' @title Choose lags of an object of class \code{timedom}
#' @param A an object of class \code{timedom}.
#' @param lags a vector which contains a set of lags. These lags must be a subset of the lags defined for timedom object A. Only those lags will be kept, the other lags are removed.
#' @return An object of class \code{timedom}.
#' @keywords time.domain
#' @export
timedom.trunc = function(A, lags){
  if (!is.timedom(A))
    stop ("A must be an object of class timedom")
   suboperators=A$operators[,,A$lags %in% lags,drop=FALSE]
   #drop=False guarantees that filter of length one is not converted from array to matrix
   timedom(suboperators,lags=as.vector(intersect(lags,A$lags)))
  }

