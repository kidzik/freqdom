#' Checks if an object belongs to the class \code{\link{timedom}}.
#'
#' @title Checks if an object belongs to the class timedom
#' @param X some object
#' @return \code{TRUE} if \code{X} is of type \code{\link{timedom}}, \code{FALSE} otherwise
#' @seealso \code{\link{freqdom}}, \code{\link{timedom}}, \code{\link{is.freqdom}}
#' @export
#' @keywords classes
is.timedom = function (X){
  inherits(X,'timedom')
}
