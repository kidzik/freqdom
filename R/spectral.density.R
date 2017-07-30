#' Estimates the spectral density and cross spectral density of vector time series.
#'
#' Let \eqn{[X_1,\ldots, X_T]^\prime} be a \eqn{T\times d_1} matrix and \eqn{[Y_1,\ldots, Y_T]^\prime} be a \eqn{T\times d_2} matrix. We stack the vectors and assume that \eqn{(X_t^\prime,Y_t^\prime)^\prime} is a stationary multivariate time series of dimension \eqn{d_1+d_2}. The cross-spectral density between the two time series \eqn{(X_t)} and \eqn{(Y_t)} is defined as
#' \deqn{
#'   \sum_{h\in\mathbf{Z}} \mathrm{Cov}(X_h,Y_0) e^{-ih\omega}.
#' }
#' The function \code{spectral.density} determines the empirical cross-spectral density between the two time series \eqn{(X_t)} and \eqn{(Y_t)}. The estimator is of form
#' \deqn{
#'   \widehat{\mathcal{F}}^{XY}(\omega)=\sum_{|h|\leq q} w(|k|/q)\widehat{C}^{XY}(h)e^{-ih\omega},
#' }
#' with \eqn{\widehat{C}^{XY}(h)} defined in \code{cov.structure} Here \eqn{w} is a kernel of the specified type and \eqn{q} is the window size.  By default the Bartlett kernel \eqn{w(x)=1-|x|} is used.
#'
#' See, e.g., Chapter 10 and 11 in Brockwell and Davis (1991) for details.
#'
#' @title Compute empirical spectral density
#' @param X a vector or a vector time series given in matrix form. Each row corresponds to a timepoint.
#' @param Y a vector or vector time series given in matrix form. Each row corresponds to a timepoint.
#' @param freq a vector containing frequencies in \eqn{[-\pi, \pi]} on which the spectral density should be evaluated.
#' @param q window size for the kernel estimator, i.e. a positive integer.
#' @param weights kernel used in the spectral smoothing. By default the Bartlett kernel is chosen.
#' @return Returns an object of class \code{\link{freqdom}}. The list is containing the following components:
#' * \code{operators} \eqn{\quad} an array. The \eqn{k}-th matrix in this array corresponds to the spectral density matrix evaluated at the \eqn{k}-th frequency listed in \code{freq}.
#' * \code{freq} \eqn{\quad} returns argument vector \code{freq}.
#' @references Peter J. Brockwell and Richard A. Davis
#' \emph{Time Series: Theory and Methods}
#' Springer Series in Statistics, 2009
# @noRd
#' @export
#' @keywords dpca
spectral.density = function(X, Y = X,
                            freq = (-1000:1000/1000)*pi,
                            q = max(1,floor(dim(X)[1]^(1/3))),
                            weights = c('Bartlett', 'trunc', 'Tukey', 'Parzen', 'Bohman', 'Daniell', 'ParzenCogburnDavis')){
  if (length(weights)>1)
    weights = "Bartlett"
  if (is.vector(X)) X=as.matrix(X)
  if (is.vector(Y)) Y=as.matrix(Y)
  if (!is.matrix(X) || !is.matrix(Y))
    stop("X and Y must be matrices")
  if (dim(X)[1] != dim(Y)[1])
    stop("Number of observations must be equal")
  if (!(q > 0))
    stop("q must be a positive integer")

  thetas = freq

  nbasisX = dim(X)[2]
  nbasisY = dim(Y)[2]
  n = dim(X)[1]
  Ch = cov.structure(X,Y,-q:q)

  if (weights=="Bartlett"){
    wfunc = weights.Bartlett
  }else if (weights=="trunc"){
    wfunc = weights.trunc
  }else if (weights=="Tukey"){
    wfunc = weights.Tukey
  }else if (weights=="Parzen"){
    wfunc = weights.Parzen
  }else if (weights=="Bohman"){
    wfunc = weights.Bohman
  }else if (weights=="Daniell"){
    wfunc = weights.Daniell
  }else if (weights=="ParzenCogburnDavis"){
    wfunc = weights.ParzenCogburnDavis
  }else
    stop(paste(weights, "is not a valid weight function"))

  w = wfunc(-q:q/q)

  for (i in 1:dim(Ch$operators)[3])
    Ch$operators[,,i] = w[i] * Ch$operators[,,i]

  fourier.transform(Ch, freq=thetas)
}
