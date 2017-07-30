#' For a given spectral density matrix dynamic principal component filter sequences are computed.
#'
#' Dynamic principal components are linear filters \eqn{(\phi_{\ell k}\colon k\in \mathbf{Z})},
#' \eqn{1 \leq \ell \leq d}. They are defined as the Fourier coefficients of the dynamic eigenvector
#' \eqn{\varphi_\ell(\omega)} of a spectral density matrix \eqn{\mathcal{F}_\omega}:
#' \deqn{
#'   \phi_{\ell k}:=\frac{1}{2\pi}\int_{-\pi}^\pi \varphi_\ell(\omega) \exp(-ik\omega) d\omega.
#' }
#' The index \eqn{\ell} is referring to the \eqn{\ell}-th #'largest dynamic eigenvalue. Since the \eqn{\phi_{\ell k}} are
#' real, we have \deqn{
#' \phi_{\ell k}^\prime=\phi_{\ell k}^*=\frac{1}{2\pi}\int_{-\pi}^\pi \varphi_\ell^* \exp(ik\omega)d\omega.
#' }
#' For a given
#' spectral density (provided as on object of class \code{freqdom}) the function
#' \code{dpca.filters()} computes \eqn{(\phi_{\ell k})} for \eqn{|k| \leq} \code{q}  and \eqn{1 \leq \ell \leq} \code{Ndpc}.
#'
#' For more details we refer to Chapter 9 in Brillinger (2001), Chapter 7.8 in Shumway and
#' Stoffer (2006) and to Hormann et al. (2015).
#'
#' @title Compute DPCA filter coefficients
#' @param F \eqn{(d\times d)} spectral density matrix, provided as an object of class \code{freqdom}.
#' @param Ndpc an integer \eqn{\in\{1,\ldots, d\}}. It is the number of dynamic principal
#' components to be computed. By default it is set equal to \eqn{d}.
#' @param q a non-negative integer. DPCA filter coefficients at lags \eqn{|h|\leq} \code{q} will be computed.
#' @return An object of class \code{timedom}.  The list has the following components:
#' * \code{operators} \eqn{\quad} an array. Each matrix in this array has dimension \code{Ndpc} \eqn{\times d} and is
#' assigned to a certain lag. For a given lag \eqn{k}, the rows of the matrix correpsond to
#' \eqn{\phi_{\ell k}}.
#' * \code{lags} \eqn{\quad} a vector with the lags of the filter coefficients.
#' @references Hormann, S., Kidzinski, L., and Hallin, M.
#' \emph{Dynamic functional principal components.} Journal of the Royal
#' Statistical Society: Series B (Statistical Methodology) 77.2 (2015): 319-348.
#' @references Brillinger, D.
#' \emph{Time Series} (2001), SIAM, San Francisco.
#' @references Shumway, R.H., and Stoffer, D.S.
#' \emph{Time Series Analysis and Its Applications} (2006), Springer, New York.
#' @seealso \code{\link{dpca.var}}, \code{\link{dpca.scores}}, \code{\link{dpca.KLexpansion}}
#' @keywords DPCA
#' @export
dpca.filters = function(F, Ndpc = dim(F$operators)[1], q = 30){
  if (!is.freqdom(F))
    stop("F must be an object of class freqdom")
  if (q <= 0)
    stop("q must be positive")

  d=dim(F$operators)[1]
  lags_full = -q:q
  E = freqdom.eigen(F)

  XI = array(0, c(d, d, length(lags_full)))

  fd.E = freqdom(E$vectors, freq = F$freq)

  A = fourier.inverse(fd.E, lags = lags_full)

  # determine the number of lags from the threshold
  #if (!is.null(thresh)){
  #  nms = timedom.norms(A)$norms
  #  r = c(nms[q+1],nms[q:1] + nms[(q+2):(2*q+1)])
  #  rsum = sum(r)
  #   w = which(cumsum(r / rsum) > thresh)
  #  if (length(w) == 0)
  #    lags = lags_full
  #  else
  #    lags = (-(w[1]-1)):(w[1]-1)
  #  A = timedom.trunc(A, lags = lags)
  #}
  if (Ndpc < dim(A$operators)[1]){
    newdim = dim(A$operators)
    newdim[1] = Ndpc
    A$operators = array(A$operators[1:Ndpc,,],newdim)
  }
  timedom(A$operators,lags=lags_full)
}
