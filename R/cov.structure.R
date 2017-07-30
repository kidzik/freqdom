#' This function computes the empirical cross-covariance of two stationary multivariate time series.
#' If only one time series is provided it determines the empirical autocovariance function.
#'
#' Let \eqn{[X_1,\ldots, X_T]^\prime} be a \eqn{T\times d_1} matrix and
#' \eqn{[Y_1,\ldots, Y_T]^\prime} be a \eqn{T\times d_2} matrix. We stack the vectors and
#' assume that \eqn{(X_t^\prime,Y_t^\prime)^\prime} is a stationary multivariate time series
#' of dimension \eqn{d_1+d_2}. This function determines empirical lagged covariances between
#' the series \eqn{(X_t)} and \eqn{(Y_t)}. More precisely it determines
#' \eqn{\widehat{C}^{XY}(h)} for \eqn{h\in} lags,
#' where \eqn{\widehat{C}^{XY}(h)} is the empirical version of \eqn{\mathrm{Cov}(X_h,Y_0)}.
#' For a sample of size \eqn{T} we set \eqn{\hat\mu^X=\frac{1}{T}\sum_{t=1}^T X_t$ and $\hat\mu^Y=\frac{1}{T}\sum_{t=1}^T Y_t} and
#' \deqn{\hat C^{XY}(h) = \frac{1}{T}\sum_{t=1}^{T-h} (X_{t+h} - \hat\mu^X)(Y_{t} - \hat\mu^Y)'}
#' and for \eqn{h < 0}
#' \deqn{\hat C^{XY}(h) = \frac{1}{T}\sum_{t=|h|+1}^{T} (X_{t+h} - \hat\mu^X)(Y_{t} - \hat\mu^Y)'.}
#'
#' @title Estimate cross-covariances of two stationary multivariate time series
#'
#' @param X vector or matrix. If matrix, then each row corresponds to a timepoint of a vector time series.
#' @param Y vector or matrix. If matrix, then each row corresponds to a timepoint of a vector time series.
#' @param lags an integer-valued vector \eqn{(\ell_1,\ldots, \ell_K)} containing the lags for which covariances are calculated.
#' @return An object of class \code{\link{timedom}}. The list contains
#' * \code{operators} \eqn{\quad} an array. Element \code{[,,k]} contains the covariance matrix related to lag \eqn{\ell_k}.
#' * \code{lags} \eqn{\quad} returns the lags vector from the arguments.
# @references Peter J. Brockwell and Richard A. Davis
# \emph{Time Series: Theory and Methods}
# Springer Series in Statistics, 2009
#' @keywords time.domain
#' @export

cov.structure = function(X,Y=X,lags=0){
  # if no Y compute autocovariance
	if (is.null(Y))
		Y = X
	if (is.vector(X)) X = as.matrix(X)
	if (is.vector(Y)) Y = as.matrix(Y)
	lrange = as.integer(max(abs(lags)))

	if (!is.integer(lrange))
	  stop("lags must be an integer")
	if (!is.matrix(X) || !is.matrix(Y))
	  stop("X and Y must be matrices")
	if (dim(X)[1] != dim(Y)[1])
	  stop("Number of observations must be equal")

	nbasisX = dim(X)[2]
	nbasisY = dim(Y)[2]
	n = dim(X)[1]

	Ch = array(0,c(nbasisX,nbasisY,2*lrange+1))

	for (h in (-lrange):lrange)
		Ch[,,h+lrange+1] = lagged.cov(X,Y,h)
	Ch

  A = timedom(Ch,-lrange:lrange)
  timedom.trunc(A, lags)
}
