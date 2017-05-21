#' For given scores \code{Y} and dynamic principal components \code{XI}
#' retrive a series from which scores \code{Y} were calculated.
#' This procedure should be seen as the inverse of \code{\link{dpca.scores}}.
#'
#' @title Retrieve a process from given scores
#' @param Y scores process
#' @param XI principal components series
#' @return Retrived process X
#' @seealso \code{\link{dpca.scores}}, \code{\link{dprcomp}}
#' @references Hormann Siegfried, Kidzinski Lukasz and Hallin Marc.
#' \emph{Dynamic functional principal components.} Journal of the Royal
#' Statistical Society: Series B (Statistical Methodology) 77.2 (2015): 319-348.
# @export
dpca.inverse = function(Y,XI){
  t(rev(XI)) %c% Y
}
