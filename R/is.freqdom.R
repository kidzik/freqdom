inClass = function(X,cls){
  !is.null(oldClass(X)) && oldClass(X) == cls
}

#' Checks if an object belongs to the class (\code{\link{freqdom}}).
#'
#' @title Checks if an object belongs to the class freqdom
#' @param X some object
#' @return \code{TRUE} if \code{X} is of type \code{\link{freqdom}}, \code{FALSE} otherwise
#' @seealso \code{\link{freqdom}}, \code{\link{timedom}}, \code{\link{is.timedom}}
#' @export 
is.freqdom = function (X){
  inClass(X,'freqdom')
}

is.positiveint = function (n){
  is.numeric(n) && n > 0
}
