#' * For a given \code{freqdom} object \code{x}, the function \code{rev} reverts the order of the evaluation grid.
#' 
#' @name rev
#' @title Invert order of lags or grid parameters of a \code{\link{timedom}} or \code{\link{freqdom}} object, respectively
#' @param x an object of class \code{freqdom} or \code{timedom}.
#' @return Returns object of same class as \code{x}.
# @describeIn rev Reverts order of frequency grid in an object of class \code{\link{freqdom}}
#' @keywords internal
#' @export
rev.freqdom = function(x){
  order=length(x$freq):1
  x$operators = x$operators[,,order,drop=FALSE]
  dimnames(x$operators)[[3]]<-paste("frequency", x$freq)
  x
}

#' * For a given \code{timedom} object \code{x}, the function \code{rev} reverts the order of lags.
#'
#' @title rev Reverts order of lags in an object of class \code{\link{timedom}}
#' @keywords internal
#' @export
rev.timedom = function(x){
  order=length(x$lags):1
  x$operators = x$operators[,,order,drop=FALSE]
  dimnames(x$operators)[[3]]<-paste("lag", x$lags)
  x
}