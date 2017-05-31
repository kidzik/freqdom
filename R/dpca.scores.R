#' Computes the dynamic principal component scores of a vector time series.
#' 
#' The \eqn{\ell}-th dynamic principal components score sequence is defined by
#' \deqn{
#'   Y_{\ell t}:=\sum_{k\in\mathbf{Z}}  \phi_{\ell k}^\prime X_{t-k},\quad 1\leq \ell\leq d,
#' }
#' where \eqn{\phi_{\ell k}} are the dynamic PC filters as explained in  \code{dpca.filters}. For the sample version the sum extends over the range of lags for which the \eqn{\phi_{\ell k}} are defined. The actual operation carried out is \code{filter.process(X, A = dpcs)}.
#' 
#' We for more details we refer to Chapter 9 in Brillinger (2001), Chapter 7.8 in Shumway and Stoffer (2006) and to H"ormann et al. (2015).
#'
#' @title Compute scores of dynamic principal components
#' @param X a vector time series given as a \eqn{(T\times d)}-matix. Each row corresponds to a timepoint.
#' @param dpcs an object of class \code{timedom}, representing the dpca filters obtained from the sample X. If \code{dpsc = NULL}, then \code{dpcs =
#' dpca.filter(spectral.density(X))} is used.
#' @return A \eqn{(T\times Ndpc)}-matix with \code{Ndpc = dim(dpcs\eqn{\$}operators)[1]}. The \eqn{\ell}-th column contains the \eqn{\ell}-th dynamic principal component score sequence.
#' @seealso \code{\link{dpca.inverse}}, \code{\link{dprcomp}}
#' @references Hormann Siegfried, Kidzinski Lukasz and Hallin Marc.
#' \emph{Dynamic functional principal components.} Journal of the Royal
#' Statistical Society: Series B (Statistical Methodology) 77.2 (2015): 319-348.
#' @references Brillinger, D.
#' \emph{Time Series} (2001), SIAM, San Francisco.
#' @references Shumway, R.H., and Stoffer, D.S.
#' \emph{Time Series Analysis and Its Applications} (2006), Springer, New York.
#' @export
dpca.scores = function(X,dpcs){
  dpcs %c% X
}

