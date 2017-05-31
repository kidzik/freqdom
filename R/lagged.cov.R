#' For given two multivariate stationary time series \eqn{X} and \eqn{Y} function \code{lagged.cov}
#' estimates the lagged covariance at lag \eqn{h}.
#' 
#' Given a stationary multivariate time series (\eqn{X_t}',\eqn{Y_t}')' of dimentian \eqn{p_1 + p_2}
#' and a lag \eqn{h} this function determines empirical covariance between the series
#' \eqn{X_{t + h}} and \eqn{Y_t}. More precisely, it returns the empirical version of
#' \eqn{Cov(X_{h} Y_0)}. 
#' For a sample of size \eqn{T} we set \eqn{\mu^X = \frac{1}{T}\sum_{t=1}^T X_t} and
#' \eqn{\mu^Y = \frac{1}{T}\sum_{t=1}^T Y_t}. Then for \eqn{h \geq 0} we define
#' \deqn{\hat C^{XY}(h) = \frac{1}{T}\sum_{t=1}^{T-h} (X_{t+h} - \hat\mu^X)(Y_{t} - \hat\mu^Y)'}
#' and for \eqn{h < 0}
#' \deqn{\hat C^{XY}(h) = \frac{1}{T}\sum_{t=|h|+1}^{T} (X_{t+h} - \hat\mu^X)(Y_{t} - \hat\mu^Y)'.}
#' 
#' If \eqn{Y_t} is not given, we assume \eqn{Y_t := X_t} and the function will return a lagged autocovariance of \eqn{X_t}.
#' 
#' Use \code{\link{cov.structure}} to estimate the entire covariance structure on a 
#' predefined set of lags.
#'
#' @title Estimates the lagged covariance of two multivariate stationary time series
#'
#' @param X a multivariate time series represented as a \eqn{T \times p_1} matrix,
#'  where \eqn{T} is the number of observations and \eqn{p_1} is the number of covariates.
#' @param Y a multivariate time series represented as a \eqn{T \times p_2} matrix,
#'  where \eqn{T} is the number of observations and \eqn{p_2} is the number of covariates.
#'  If \eqn{Y = NULL} then \eqn{Y := X} is used and an autocovariance of \eqn{X} are computed.
#' @param lag a requested lag to evaluate the lagged covariance.
#' @return Function returns a matrix in \eqn{\mathbf{R}^{p_1 \times p_2}} corresponding to the lagged covariance \eqn{Cov(X_{h} Y_0)}, where \eqn{h} is the requested lag.
#' @seealso \code{\link{cov.structure}}
#' @references Peter J. Brockwell and Richard A. Davis
#' \emph{Time Series: Theory and Methods}
#' Springer Series in Statistics, 2009
#' @noRd
# @export
#' @examples
#' X = rar(100)
#' Y = rar(100)
#' lagged.cov(X,Y)
lagged.cov = function(X,Y=NULL,lag=0){
  if (base::is.null(Y))
		Y = X

  if (base::dim(X)[1] != base::dim(Y)[1])
    stop("Number of observations must be equal")
  if (!base::is.matrix(X) || !base::is.matrix(Y))
    stop("X and Y must be matrices")
  
	n = base::dim(X)[1]
	h = base::abs(lag)
	
  if (n - 1 <= h)
	  base::stop(base::paste("Too little observations to compute lagged covariance with lag",h))
	
  M = base::t(X[1:(n-h),]) %*% (Y[1:(n-h)+h,])/(n)
	if (lag < 0){
	  M = base::t(Y[1:(n-h),]) %*% (X[1:(n-h)+h,])/(n)
	  M = base::t(M)
	}
	M
}

