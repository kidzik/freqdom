#' For a given frequency-domain operator \code{S} (\code{\link{freqdom}}) the function \code{freqdom.transpose} computes the transpose of \eqn{S_\theta'} at
#' each frequency from the evaluation grid.
#' 
#' Let \eqn{S = \{ S_\theta : \theta \in G \}}, where \eqn{G} is some finite grid
#' of frequencies in \eqn{[-\pi,\pi]} and \eqn{S_\theta \in \mathbf{C}^{p \times p}}.
#' At each frequency \eqn{\theta \in G} function \code{freqdom.transpose} transposes
#'
#' Resulting object is defined as
#' \deqn{S' = \{ S_\theta': \theta \in G \}.}
#' 
#' @title Compute a transpose of a given frequency-domain operator at each frequency
#' @param S a frequency-domain filter of type \code{\link{freqdom}}, i.e. a set of linear operators \eqn{S_k \in \mathbf{R}^{p_1 \times p_2}}
#' on some discreet grid defined of \eqn{[-\pi,\pi]}. 
#' @return Function returns a frequency domain object (\code{\link{freqdom}}) of dimensions \eqn{L \times p_2 \times p_1}, where \eqn{L} is the size of the grid.
#' The elements of the object correspond to \eqn{S_\theta'} as defined above.
#' @seealso \code{\link{freqdom.product}}, \code{\link{freqdom.ratio}}
#' @export
freqdom.transpose = function(x){
  lags = freqdom.lags(x)
  for (i in 1:length(lags))
    x$operators[i,,] = t(x$operators[i,,])
  x
}

#' @export t freqdom
t.freqdom = freqdom.transpose

#' @export t timedom
t.timedom = freqdom.transpose
