#' For a given frequency-domain operators \code{F} and \code{G} (\code{\link{freqdom}}) the function \code{freqdom.ratio} computes \eqn{F G_\theta^{-1}} at
#' each frequency from the evaluation grid. This is a simple convenience function returning the \code{freqdom.product} of \eqn{F} and \eqn{H}, where
#' \eqn{H} is the \code{freqdom.inverse} of \eqn{G}.
#' 
#' @title Compute an inverse of a given frequency-domain operator at each frequency
#' @param F frequency-domain filter of type \code{\link{freqdom}}, i.e. a set of linear operators \eqn{F_k \in \mathbf{R}^{p_1 \times p_2}}
#' defined on discreet grid \eqn{S \subset [-\pi,\pi]}. 
#' @param G frequency-domain filter of type \code{\link{freqdom}}, i.e. a set of linear operators \eqn{G_k \in \mathbf{R}^{q_1 \times q_2}}
#' defined on a discreet grid \eqn{S \subset [-\pi,\pi]}. 
#' @param n number of observations used for estimation - the precision of inversion is calculated using this parameter. Either \code{n} or \code{K} must be given. 
#' @param K how many directions should be inverted (as in \code{\link{pseudoinverse}}). Either \code{n} or \code{K} must be given. 
#' @return Function returns a frequency domain object (\code{\link{freqdom}}) of the same dimensions as \eqn{F}.
#' The elements of the object correspond to \eqn{F_\theta G_\theta^{-1}} as defined above.
#' @references Hormann Siegfried, Kidzinski Lukasz and Kokoszka Piotr.
#' \emph{Estimation in functional lagged regression.} 
#' Journal of Time Series Analysis 36.4 (2015): 541-561.
#' @noRd
# @export
#' @seealso \code{\link{freqdom.product}, \link{reg.dim.est}}, \code{\link{freqdom.product}}, \code{\link{freqdom.inverse}}
#' @examples
#' n = 100
#' X = rar(n)
#' Y = rar(n)
#' SYX = spectral.density(X,Y)
#' SXX = spectral.density(X)
#' R = freqdom.ratio(SYX,SXX,n)
# freqdom.ratio = function(F,G,n=NULL,K=NULL){
#   Ginv = freqdom.inverse(G,n=n,K=K)
#   F %*% Ginv
# }