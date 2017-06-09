#' * For a given frequency-domain operator \code{x} (\code{\link{freqdom}}), the function \code{rev} inverts the evaluation grid.
#' 
#' Let \eqn{S = \{ S_k : k \in G \}} be an ordered set of operators, where \eqn{G} is a finite grid
#' of frequencies in \eqn{G \subset [-\pi,\pi]} or set of lags \eqn{G \subset \mathbf{Z}}.
#' Let \eqn{k_1, k_2, ..., k_L} be a sorted sequence of elements of \eqn{G}
#' Function \code{rev} inverts the order of this sequence, i.e. it returns an ordered set of operators
#' \deqn{F = \{ F_{\theta_i} : \theta_{L - i + 1} \in G \}.}
#' 
#' @name rev
#' @title Invert the lags or the grid parameters of a \code{\link{freqdom}} or a \code{\link{timedom}} object
#' @param x a frequency-domain filter of type \code{\link{freqdom}} or a time-domain filter of type \code{\link{timedom}}, i.e. a set of linear operators \eqn{S_k \in \mathbf{R}^{p_1 \times p_2}}
#' defined on \eqn{k \in G \subset \mathbf{Z}.}
#' @return Function returns a frequency-domain object (\code{\link{freqdom}}) or a time-domain object (\code{\link{time}}) of the same dimensions,
#' with the same operators but with the inverted order of operators.
#' @export
# @describeIn rev Reverts time in a \code{\link{timedom}} or a \code{\link{freqdom}} object
# @describeIn rev Reverts time in a \code{\link{freqdom}} object
#' @export
rev.freqdom = function(x){
  order=length(x$freq):1
  x$operators = x$operators[,,order,drop=FALSE]
  dimnames(x$operators)[[3]]<-paste("frequency", x$freq)
  x
}

#' * For a given time-domain operator \code{x} (\code{\link{timedom}}), the function \code{rev} inverts time.
#' 
#' @describeIn rev Reverts time in a \code{\link{timedom}} object
#' @export
rev.timedom = function(x){
  order=length(x$lags):1
  x$operators = x$operators[,,order,drop=FALSE]
  dimnames(x$operators)[[3]]<-paste("lag", x$lags)
  x
}