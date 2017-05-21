inClass = function(X,cls){
  !is.null(oldClass(X)) && oldClass(X) == cls
}

#' Checks if a given object is of frequency domain type (\code{\link{freqdom}}).
#'
#' @title Checks if a given object is of frequency domain type
#' @param X any object
#' @return \code{TRUE} if \code{X} is of type \code{\link{freqdom}}, \code{FALSE} otherwise
#' @seealso \code{\link{freqdom}}, \code{\link{timedom}}, \code{\link{is.timedom}}
#' @export 
is.freqdom = function (X){
  inClass(X,'freqdom')
}

is.positiveint = function (n){
  is.numeric(n) && n > 0
}
