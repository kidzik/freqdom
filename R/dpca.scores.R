#' Compute scores of a given series \code{X} using dynamic principal components \code{XI}.
#' Procedure can be inverted using \code{\link{dpca.inverse}}
#'
#' @title Compute scores of dynamic principal components
#' @param X series to 'project'
#' @param XI principal components series
#' @return Matrix of scores of \code{X}
#' @seealso \code{\link{dpca.inverse}}, \code{\link{dprcomp}}
#' @references Hormann Siegfried, Kidzinski Lukasz and Hallin Marc.
#' \emph{Dynamic functional principal components.} Journal of the Royal
#' Statistical Society: Series B (Statistical Methodology) 77.2 (2015): 319-348.
# @export
dpca.scores = function(X,XI){
  XI %c% X
}

