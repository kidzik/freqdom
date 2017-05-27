#' Estimate the covariance structure of two stationary multivariate time series
#' on a given set of lags \eqn{[-q,q]}.
#' 
#' Given a stationary multivariate time series (\eqn{X_t}',\eqn{Y_t}')' of dimentian \eqn{p_1 + p_2}
#' this function determines empirical lagged covariances between the series \eqn{X_t} and \eqn{Y_t}. More
#' precisely, it returns 
#' \deqn{(\hat C^{XY}(h): h \in \{ -q,...,0,...,q\})}
#' where \eqn{\hat C^{XY}(h)} is the empirical version of \eqn{Cov(X_{h} Y_0)}. 
#' For a sample of size \eqn{T} we set \eqn{\mu^X = \frac{1}{T}\sum_{t=1}^T X_t} and
#' \eqn{\mu^Y = \frac{1}{T}\sum_{t=1}^T Y_t}. Then for \eqn{h \geq 0} we define
#' \deqn{\hat C^{XY}(h) = \frac{1}{T}\sum_{t=1}^{T-h} (X_{t+h} - \hat\mu^X)(Y_{t} - \hat\mu^Y)'}
#' and for \eqn{h < 0}
#' \deqn{\hat C^{XY}(h) = \frac{1}{T}\sum_{t=|h|+1}^{T} (X_{t+h} - \hat\mu^X)(Y_{t} - \hat\mu^Y)'.}
#' 
#' If \eqn{Y_t} is not given, we assume \eqn{Y_t := X_t} and the function will return autocovariances of \eqn{X_t}.
#' 
#' Use \code{\link{lagged.cov}} to estimate a particular lag.
#'
#' @title Estimate cross-covariances of two stationary multivariate time series
#'
#' @param X a multivariate time series represented as a \eqn{T \times p_1} matrix,
#'  where \eqn{T} is the number of observations and \eqn{p_1} is the number of covariates.
#' @param Y a multivariate time series represented as a \eqn{T \times p_2} matrix,
#'  where \eqn{T} is the number of observations and \eqn{p_2} is the number of covariates.
#'  If \eqn{Y = NULL} then \eqn{Y := X} is used and autocovariances of \eqn{X} are computed.
#' @param q covariances for lags \eqn{|h| \leq q} 
#' \code{q} must be a positive integer.
#' @return Function returns a time domain object (\code{\link{timedom}}) of dimensions \eqn{(2q + 1) \times p_1 \times p_2}
#' representing the estimated covariance structure on lags \eqn{\{ -q,...,0,...,q\}}. 
#' @seealso \code{\link{lagged.cov}}
#' @references Peter J. Brockwell and Richard A. Davis
#' \emph{Time Series: Theory and Methods}
#' Springer Series in Statistics, 2009
#' @export
#' @examples
#' X = rar(100)
#' Y = rar(100)
#' cov.structure(X,Y)
cov.structure = function(X,Y=X,lags=0){
  # if no Y compute autocovariance
	if (is.null(Y))
		Y = X
	
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

	Ch = array(0,c(2*lrange+1,nbasisX,nbasisY))
	
	for (h in (-lrange):lrange)
		Ch[h+lrange+1,,] = lagged.cov(X,Y,h)
	Ch
  
  A = timedom(Ch,-lrange:lrange)
  timedom.trunc(A, lags)
}
