#' Computes the dynamic Karhunen-Loeve expansion of a vector time series up to a given order.
#' 
#' We obtain the dynamic Karhnunen-Loeve expansion of order \eqn{L}, \eqn{1\leq L\leq d}. It is defined as
#' \deqn{
#'   \sum_{\ell=1}^L\sum_{k\in\mathbf{Z}} Y_{\ell, t+k} \phi_{\ell k},
#' }
#' where \eqn{\phi_{\ell k}} are the dynamic PC filters as explained in  \code{\link{dpca.filters}}. For the sample version the sum in \eqn{k} extends over the range of lags for which the \eqn{\phi_{\ell k}} are defined. The actual operation carried out is the command
#' \code{filter.process(dpca.scores(X, dpcs),t(rev(dpsc)))}. The function \code{\link{rev}} reverts time.
#' 
#' For more details we refer to Chapter 9 in Brillinger (2001), Chapter 7.8 in Shumway and Stoffer (2006)
#' and to Hormann et al. (2015).
#'
#' @title Retrieve a process from given scores
#' @param X a vector time series given as a \eqn{(T\times d)}-matix. Each row corresponds to a timepoint.
#' @param dpcs an object of class \code{timedom}, representing the dpca filters obtained from the sample X. If \code{dpsc = NULL},
#' then \code{dpcs = dpca.filter(spectral.density(X))} is used.
#' @return A \eqn{(T\times d)}-matix. The \eqn{\ell}-th column contains the \eqn{\ell}-th data point.
#' @references Hormann Siegfried, Kidzinski Lukasz and Hallin Marc.
#' \emph{Dynamic functional principal components.} Journal of the Royal
#' Statistical Society: Series B (Statistical Methodology) 77.2 (2015): 319-348.
#' @references Brillinger, D.
#' \emph{Time Series} (2001), SIAM, San Francisco.
#' @references Shumway, R.H., and Stoffer, D.S.
#' \emph{Time Series Analysis and Its Applications} (2006), Springer, New York.
#' @seealso \code{\link{dpca.filters}}, \code{\link{filter.process}}, \code{\link{dpca.scores}}, \code{\link{rev}}
#' @export
dpca.KLexpansion = function(X,dpcs){
  t(rev(dpcs)) %c% X
}
