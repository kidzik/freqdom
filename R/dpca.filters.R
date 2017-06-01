#' For a given spectral density matrix dynamic principal component filter sequences are computed.
#' 
#' Dynamic principal components are linear filters \eqn{(\phi_{\ell k}\colon k\in \mathbf{Z})}, \eqn{1 \leq \ell \leq d}. They are defined as the Fourier coefficients of the dynamic eigenvector \eqn{\varphi_\ell(\omega)} of a spectral density matrix \eqn{\mathcal{F}_\omega}:
#' \deqn{
#'   \phi_{\ell k}:=\frac{1}{2\pi}\int_{-\pi}^\pi \varphi_\ell(\omega) \exp(-ik\omega) d\omega.
#' }
#' The index \eqn{\ell} is referring to the \eqn{\ell}-th largest dynamic eigenvalue. For a given spectral density (provided as on object of class \code{freqdom}) the function \code{dpca.filters()} computes \eqn{(\phi_{\ell k}\colon |k| \leq q)}. Filters will be computed for \eqn{1 \leq \ell \leq} Ndpc.
#' 
#' Note that \eqn{\sum_{k\in\mathbf{Z}}\|\phi_{\ell k}\|^2=1}. If q is not provided, then q will be chosen as the smallest integer such that \eqn{\sum_{|k| \leq q}\|\phi_{\ell k}\|^2 \geq} thresh.
#' 
#' For more details we refer to Chapter 9 in Brillinger (2001), Chapter 7.8 in Shumway and Stoffer (2006) and to Hormann et al. (2015).
#'
#' @title Compute DPCA filter coefficients
#' @param F \eqn{(d\times d)} spectral density matrix, provided as an object of class \code{freqdom}.
#' @param Ndpc an integer \eqn{\in\{1,\ldots, d\}}. It is the number of dynamic principal components to be computed. By default it is set equal to \eqn{d}.
#' @param q a non-negative integer. DPCA filter coefficients at lags \eqn{|h|\leq} q will be computed.
#' @param thresh a numerical number in \eqn{(0,1)}. Can be provided alternatively to q. The closer to 1, the more filter coefficients will be computed.
#' @return An object of class \code{timedom}.  The list has the following components:
#' * \code{operators} an array. Each matrix in this array has dimension Ndpc \eqn{\times d} and is assigned to a certain lag. For a given lag \eqn{k}, the rows of the matrix correpsond to \eqn{\phi_{\ell k}^\prime}.
#' * \code{lags} a vector with the lags of the filter coefficients.
#' @references Hormann Siegfried, Kidzinski Lukasz and Hallin Marc.
#' \emph{Dynamic functional principal components.} Journal of the Royal
#' Statistical Society: Series B (Statistical Methodology) 77.2 (2015): 319-348.
#' @references Brillinger, D.
#' \emph{Time Series} (2001), SIAM, San Francisco.
#' @references Shumway, R.H., and Stoffer, D.S.
#' \emph{Time Series Analysis and Its Applications} (2006), Springer, New York.
#' @seealso \code{\link{dpca.inverse}}, \code{\link{dpca.scores}}, \code{\link{dpca.KLexpansion}}
#' @export
dpca.filters = function(F, L = 10){
  if (!is.freqdom(F))
    stop("F must be a freqdom operator")
  if (L < 0)
    stop("L must be a non-negative integer")

  lags = -L:L
  E = freqdom.eigen(F)
  
	d = dim(E$vectors)[2]

  XI = array(0,c(length(lags),d,d))

	for (component in 1:d)
    XI[,component,] = t(exp(-(F$freq %*% t(lags)) * 1i)) %*% E$vectors[,,component] / length(F$freq) 
  # TODO: this is just 'fourier.inverse' written in a fancy way - should be simplified

	timedom(Re(XI[length(lags):1,,]),lags)
}
