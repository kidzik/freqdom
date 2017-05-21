rev.freqdom = function(x){
  x$freq = rev(x$freq)
  x
}

#' @export
rev.timedom = function(x){
  x$lags = rev(x$lags)
  x
}

#' For a given frequency-domain operator \code{x} (\code{\link{freqdom}}), the function \code{rev} inverts the order of the operators.
#' 
#' Let \eqn{S = \{ S_\theta : \theta \in G \}}, where \eqn{G} is a finite grid
#' of frequencies in \eqn{G \subset [-\pi,\pi]}.
#' Let \eqn{\theta_1, \theta_2, ..., \theta_L} be a sorted sequence of elements of \eqn{G}
#' Function \code{rev} inverts the order of this sequence, i.e. it returns an operator
#' \deqn{F = \{ F_{\theta_i} : \theta_{L - i + 1} \in G \}.}
#' 
#' @title Invert the lags or the grid parameters of a \code{\link{freqdom}} object
#' @param x a frequency-domain filter of type \code{\link{freqdom}}, i.e. a set of linear operators \eqn{S_k \in \mathbf{R}^{p_1 \times p_2}}
#' on a discreet grid defined of \eqn{[-\pi,\pi]}. 
#' @return Function returns a frequency domain object (\code{\link{freqdom}}) of the same dimensions and
#' with the same operators but with the inverted order of operators.
#' @export
"rev"